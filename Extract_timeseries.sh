#!/bin/bash
#
# Extract_timeseries.sh
# (C) C.Vriend - AmsUMC dec 2022 (refactored)
#
# Warps standard atlases (Schaefer, optionally Brainnetome and/or the Buckner
# cerebellar atlas) to subject-specific FreeSurfer space, resamples them into
# the subject's BOLD space, and extracts ROI timeseries + ROI volumes from the
# denoised functional MRI (output of the rsfmridenoise_sarray_v3.sh step).
#
#SBATCH --begin=04:30:00
#SBATCH --job-name=extract_timeseries
#SBATCH --mem=2G
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --time=00-0:30:00
#SBATCH --nice=2000
#SBATCH --array=1-603%20
#SBATCH --output=%x_%A_%a.log

set -euo pipefail

Usage() {
    cat <<EOF

    (C) C.Vriend - AmsUMC Oct 2026
    Warps standard atlases to subject-specific FreeSurfer space and uses them
    to extract ROI timeseries + ROI volumes from the denoised, smoothed
    functional MRI for each parcel.

    Usage: sbatch Extract_timeseries.sh <freesurferdir> <fmriprepdir> [session]

    Obligatory:
      freesurferdir = full path to the FreeSurfer directory containing subject output
      fmriprepdir   = full path to the fmriprep directory containing subject output

    Optional:
      session       = session ID of the fmriprep output, e.g. ses-T0 (omit if no sessions)

    Additional options/paths (atlas list, denoise protocol, etc.) can be
    adjusted in the CONFIG section below. Atlas list entries may carry a
    "+Buckner" suffix (e.g. 400P17N+Buckner) to merge the Buckner (2011)
    cerebellar parcellation into the cortical parcellation before timeseries
    extraction.

EOF
    exit 1
}

[ $# -lt 2 ] && Usage



######################
# to add some color to the output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m'
######################

#!/bin/bash

GROUP="anw-np"

# If we're not running under the target group, re-exec ourselves
if [ "$(id -gn)" != "$GROUP" ]; then
    exec sg "$GROUP" -c "$0 $*"
fi


########################################
# CONFIG
########################################

fs_subjects_dir=${1}    # FreeSurfer directory (user-supplied)
fmriprepdir=${2}        # fmriprep directory
sess=${3:-}             # optional session ID, e.g. ses-T0

atlasdir=/appdata/users/P042819/atlas4FreeSurfer
tempdir=/scratch/users/P042819/tmp
freesurfercontainer=/appdata/containers/freesurfer/freesurfer-8.2.0-1.sif

outputspace=T1w                # generally correct, matches denoising output space
denoise_protocol=24HMP8PhysSpikeRegGS  # protocol used in the denoising step
FWHM=6                         # smoothing kernel used in the denoising step
smoothtype=susan

# Schaefer parcellations to warp + extract. Add/remove entries here only;
# everything downstream (labelconvert, bold-space warp, timeseries
# extraction) drives off this single list.
# An entry may carry a "+Buckner" suffix (currently only meaningful for the
# 17-network parcellations, e.g. 400P17N+Buckner) to append the Buckner (2011)
# cerebellar parcellation to the cortical parcellation: the Buckner label IDs
# are offset past the Schaefer node IDs and both atlases are merged into a
# single parcellation (space-T1w_atlas-400P17N+Buckner_dseg.nii.gz), for which
# ROI timeseries and volumes are then extracted. List the plain atlas (e.g.
# 400P17N) as an extra entry if you also want the cortical-only output.
SCHAEFER_ATLASES=(300P17N 400P17N 400P17N+Buckner)

# The Brainnetome atlas (BNA) is expensive to warp (mris_ca_label /
# mri_ca_label / mri_aparc2aseg take real time) and was NOT actually used
# downstream in the original script (the extraction loop only ever processed
# 300P7N/400P17N), so it ran on every subject for nothing. Set this to true
# to warp + extract BNA as well.
INCLUDE_BNA=false

# Buckner (2011) cerebellar parcellation. The atlas volume and the
# FreeSurfer-conformed MNI152 template are taken from the FreeSurfer
# distribution (FREESURFER_HOME/average/...); the atlas is warped to
# subject space with a SynthMorph registration of the MNI template to the
# subject's norm.mgz (see warp_buckner_atlas). Only the label lookup table
# lives in atlasdir; its IDs must match the FreeSurfer-distributed
# Buckner2011_17Networks volume.
buckner_lut=${atlasdir}/Buckner/Buckner_17Networks_orig.txt

# threads for mri_synthmorph; defaults to the number of CPUs requested
# from SLURM (keep this in sync with --cpus-per-task above)
nthreads=${SLURM_CPUS_PER_TASK:-4}


# atlas short-name -> Schaefer "order" string used in the fsaverage annot filename
declare -A SCHAEFER_ORDER=(
    [100P7N]="100Parcels_7Networks"
    [200P7N]="200Parcels_7Networks"
    [300P7N]="300Parcels_7Networks"
    [300P17N]="300Parcels_17Networks"
    [400P7N]="400Parcels_7Networks"
    [400P17N]="400Parcels_17Networks"
)

# atlas short-name -> label-ID prefix used to locate _orig.txt / _modified.txt
# lookup tables for labelconvert / header2timeseries.py
# NOTE: combined "BASE+Buckner" atlases need no entry here; their lookup
# tables are generated on the fly (see make_buckner_label_files).
declare -A ATLAS_LABEL_ID=(
    [100P7N]="Schaefer_100P7N"
    [200P7N]="Schaefer_200P7N"
    [300P7N]="Schaefer_300P7N"
    [300P17N]="Schaefer_300P17N"
    [400P7N]="Schaefer_400P7N"
    [400P17N]="Schaefer_400P17N"
    [aparc500]="aparc500_labels"
    [BNA]="BNA_labels"
    ["BNA+cerebellum"]="BNA+CER_labels"
)

mkdir -p "${tempdir}"

if [ -z "${sess}" ]; then
    sessionpath="/"
    sessionfile="_"
else
    sessionpath="/${sess}/"
    sessionfile="_${sess}_"
fi

########################################
# HELPER FUNCTIONS
########################################

load_modules() {
    module load mrtrix3 fsl freesurfer anaconda3
    export APPTAINER_BINDPATH="/appdata/users/${USER},/net/beegfs/users/${USER},/scratch/users/${USER}"

    # FreeSurfer's setup script (run by the module load above) sets
    # SUBJECTS_DIR=$FREESURFER_HOME/subjects as a side effect, clobbering
    # any prior value. Force it back to the user-supplied directory here,
    # AFTER all modules have loaded, or steps like ensure_fsfiles() will
    # rsync the read-only central fsaverage install into itself and fail
    # with "failed to modify permissions ... Operation not permitted".
    SUBJECTS_DIR="${fs_subjects_dir}"
    export SUBJECTS_DIR
}

select_subject() {
    cd "${fmriprepdir}"
    # NOTE: previously hardcoded to "sub-ARRIBA5016*" (leftover test/debug
    # code) which silently restricted the whole array to a handful of
    # subjects. Now lists all BIDS-formatted subject folders.
    ls -d sub-*/ | sed 's:/.*::' > subjects.txt

    subj=$(sed "${SLURM_ARRAY_TASK_ID}q;d" subjects.txt)

    echo -e "${BLUE}INITIALIZING...(wait a sec)${NC}"
    sleep $((RANDOM % 20 + 2))

    # extra check because the slurm array sometimes behaves unexpectedly.
    # NOTE: the original had a missing semicolon here ("... exit" was parsed
    # as an argument to echo, so the script never actually exited on a bad ID)
    if [[ "${subj}" != sub-* ]]; then
        echo -e "${RED}ERROR: '${subj}' is not a valid BIDS subject ID; exiting${NC}" >&2
        exit 1
    fi
}

check_denoised_image() {
    denoisedimg="${fmriprepdir}/${subj}${sessionpath}func/${subj}${sessionfile}task-rest_space-${outputspace}_rec-${smoothtype}${FWHM}mm_desc-${denoise_protocol}_bold.nii.gz"
    if [[ ! -f "${denoisedimg}" ]]; then
        echo -e "${RED}ERROR: smoothed/denoised image (protocol=${denoise_protocol}) not found for ${subj}; exiting${NC}" >&2
        exit 1
    fi
}

# Resolve the _orig.txt / _modified.txt label lookup tables for a given atlas.
# Prints "orig_path modified_path" (consumed via `read -r orig mod < <(...)`)
get_label_files() {
    local atlas=$1
    local orig modified

    # combined "BASE+Buckner" atlases use a generated lookup table (base
    # parcellation + offset Buckner labels); see make_buckner_label_files()
    if [[ "${atlas}" == *+Buckner ]]; then
        local base=${atlas%%+*}
        local base_orig base_mod
        read -r base_orig base_mod < <(get_label_files "${base}")
        local combined="${fmriprepdir}/${subj}${sessionpath}anat/${subj}${sessionfile}atlas-${atlas}_modified.txt"
        # regenerate the generated tables if they went missing
        # (e.g. partially cleaned anat dir)
        [[ -f "${combined}" ]] || make_buckner_label_files "${atlas}"
        echo "${base_orig} ${combined}"
        return 0
    fi

    local id=${ATLAS_LABEL_ID[${atlas}]:-}

    if [[ -z "${id}" ]]; then
        echo -e "${RED}ERROR: no label-ID mapping for atlas '${atlas}'${NC}" >&2
        exit 1
    fi

    # NOTE: fixed a stray "fi" before the "else" here; as pasted this was a
    # bash syntax error and the script would not parse at all.
    if [[ "${id}" == *Schaefer* ]]; then
        orig="${atlasdir}/Schaefer/${id}_orig.txt"
        modified="${atlasdir}/Schaefer/${id}_modified.txt"
    else
        orig="${atlasdir}/${atlas}/${id}_orig.txt"
        modified="${atlasdir}/${atlas}/${id}_modified.txt"
    fi
    echo "${orig} ${modified}"
}

ensure_fsfiles() {
    # get the .sif path from `module show freesurfer/<version>`, then:
    apptainer exec ${freesurfercontainer} cp -r /usr/local/freesurfer/8.2.0-1/subjects/fsaverage ${SUBJECTS_DIR}
    chmod -R ug+w "${SUBJECTS_DIR}/fsaverage"
    apptainer exec ${freesurfercontainer} cp \
    /usr/local/freesurfer/8.2.0-1/average/Yeo_JNeurophysiol11_MNI152/FSL_MNI152_FreeSurferConformed_1mm.nii.gz \
    /usr/local/freesurfer/8.2.0-1/average/Buckner_JNeurophysiol11_MNI152/Buckner2011_17Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz \
    ${SUBJECTS_DIR}
  

}


# Optional: warp the Brainnetome atlas to FreeSurfer space (cortex via
# mris_ca_label, subcortex via mri_ca_label), producing BNA+aseg.nii.gz.
# Only runs when INCLUDE_BNA=true.
warp_bna_atlas() {

    if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.nii.gz" ]]; then

        echo -e "${BLUE}warping Brainnetome atlas (BNA) to individual FreeSurfer space${NC}"

        if [[ ! -f "${SUBJECTS_DIR}/${subj}/label/lh.BN_Atlas.annot" ]] || \
        [[ ! -f "${SUBJECTS_DIR}/${subj}/label/rh.BN_Atlas.annot" ]]; then
            for hemi in lh rh; do
                mris_ca_label -seed 1234 -l "${SUBJECTS_DIR}/${subj}/label/${hemi}.cortex.label" \
                    "${subj}" "${hemi}" \
                    "${SUBJECTS_DIR}/${subj}/surf/${hemi}.sphere.reg" \
                    "${atlasdir}/BNA/${hemi}.BN_Atlas.gcs" "${SUBJECTS_DIR}/${subj}/label/${hemi}.BN_Atlas.annot"

                mris_anatomical_stats -mgz -cortex "${SUBJECTS_DIR}/${subj}/label/${hemi}.cortex.label" \
                    -f "${SUBJECTS_DIR}/${subj}/stats/${hemi}.BN_Atlas.stats" -b \
                    -a "${SUBJECTS_DIR}/${subj}/label/${hemi}.BN_Atlas.annot" \
                    -c "${atlasdir}/BNA/BNA_labels_orig_wsubcortex.txt" "${subj}" "${hemi}" white
            done
        fi

        if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz" ]]; then
            mri_ca_label -threads 2 "${SUBJECTS_DIR}/${subj}/mri/brain.mgz" \
                "${SUBJECTS_DIR}/${subj}/mri/transforms/talairach.m3z" "${atlasdir}/BNA/BN_Atlas_subcortex.gca" \
                "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz"

            mri_segstats --seg "${SUBJECTS_DIR}/${subj}/mri/BN_Atlas_subcortex.mgz" \
                --ctab "${atlasdir}/BNA/BNA_labels_orig_wsubcortex.txt" --excludeid 0 \
                --sum "${SUBJECTS_DIR}/${subj}/stats/BN_Atlas_subcortex.stats"
        fi

        if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.nii.gz" ]]; then
            mri_aparc2aseg --threads 2 --s "${subj}" --annot BN_Atlas --o "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.mgz"
            mrconvert "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.mgz" "${SUBJECTS_DIR}/${subj}/mri/BNA+aseg.nii.gz" -force
        fi
    fi

}


# Warp one Schaefer parcellation from fsaverage to subject FreeSurfer space
# (surface annot -> volumetric ${atlas}+aseg.nii.gz)
warp_schaefer_atlas() {
    local atlas=$1
    # strip an optional "+Buckner" suffix: the cortical parcellation is
    # warped as usual, the cerebellar atlas is merged in later on (see
    # build_buckner_combined_dseg)
    atlas=${atlas%%+*}
    local order=${SCHAEFER_ORDER[${atlas}]:-}

    if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz" ]]; then

        if [[ -z "${order}" ]]; then
            echo -e "${RED}ERROR: no Schaefer order mapping for atlas '${atlas}'${NC}" >&2
            exit 1
        fi
        echo -e "${BLUE}parcellation = ${atlas} (${order})${NC}"

        if [[ ! -f "${SUBJECTS_DIR}/${subj}/label/lh.${atlas}.annot" ]] || \
        [[ ! -f "${SUBJECTS_DIR}/${subj}/label/rh.${atlas}.annot" ]]; then
            for hemi in lh rh; do
                mri_surf2surf --srcsubject fsaverage --trgsubject "${subj}" --hemi "${hemi}" \
                    --sval-annot "${atlasdir}/Schaefer/fsaverage/label/${hemi}.Schaefer2018_${order}_order.annot" \
                    --tval "${SUBJECTS_DIR}/${subj}/label/${hemi}.${atlas}.annot"
            done
        fi

        if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" ]]; then
            mri_aparc2aseg --s "${subj}" --o "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" --annot "${atlas}"
        fi
        if [[ ! -f "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz" ]]; then
            mrconvert "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.mgz" "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz" -force
        fi

    fi   

    mri_convert "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz" --out_orientation RAS "${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz"
}


# Warp the Buckner (2011) cerebellar atlas from MNI152 space to individual
# FreeSurfer space using SynthMorph:
#   1. nonlinear registration of the FreeSurfer-conformed MNI152 template to
#      the subject's norm.mgz (moving = MNI, fixed = subject, so the warp
#      already maps MNI-space volumes into subject space; no inversion
#      needed later)
#   2. conversion of the SynthMorph warp to a FreeSurfer .m3z
#   3. application to the atlas with nearest-neighbour interpolation
#      (label images must not be interpolated)
#   4. restriction to the cerebellum (aseg 8 = L-Cerebellum-Cortex,
#      47 = R-Cerebellum-Cortex)
# The registration is cached as mri/transforms/mni_to_subj.m3z and reused on
# re-runs (it can also serve any other MNI-space volume). Only used by
# "*+Buckner" atlases.
warp_buckner_atlas() {
    local buckner_fs="${SUBJECTS_DIR}/${subj}/mri/Buckner2011_17Networks.mgz"
    [[ -f "${buckner_fs}" ]] && return 0

    local mni_t1="${SUBJECTS_DIR}/Yeo_JNeurophysiol11_MNI152/FSL_MNI152_FreeSurferConformed_1mm.nii.gz"
    local buckner_atlas="${SUBJECTS_DIR}/Buckner_JNeurophysiol11_MNI152/Buckner2011_17Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz"
    local xfm_dir="${SUBJECTS_DIR}/${subj}/mri/transforms"

    local f
    for f in "${mni_t1}" "${buckner_atlas}"; do
        [[ -f "${f}" ]] || { echo -e "${RED}ERROR: required file not found: ${f}${NC}" >&2; exit 1; }
    done

    echo -e "${BLUE}warping Buckner cerebellar atlas to individual FreeSurfer space${NC}"

    if [[ ! -f "${xfm_dir}/mni_to_subj.m3z" ]]; then
        echo -e "${BLUE}register MNI template to subject (SynthMorph)${NC}"

        # --- Step 1: nonlinear registration, MNI template -> subject norm.mgz ---
        # moving = mni_t1, fixed = subject. This way the warp already maps
        # MNI-space volumes (like the Buckner atlas) into subject space,
        # so no inversion is needed later.
        mri_synthmorph register -m joint \
            -o "${xfm_dir}/mni_in_subj_QC.mgz" \
            -t "${xfm_dir}/mni_to_subj_warp.mgz" \
            "${mni_t1}" \
            "${SUBJECTS_DIR}/${subj}/mri/norm.mgz" \
            -j "${nthreads}"

        # --- Step 2 (optional QC): eyeball the registration before trusting it ---
        # freeview -v ${SUBJECTS_DIR}/${subj}/mri/norm.mgz \
        #     ${xfm_dir}/mni_in_subj_QC.mgz:colormap=heat:opacity=0.5

        # --- Step 3: convert the SynthMorph warp to FreeSurfer .m3z ---
        # -g gives the source geometry (the moving image, i.e. mni_t1).
        # Confirm --inras is correct for your build: run
        #   mri_warp_convert --help
        # and check mri_info on mni_to_subj_warp.mgz (RAS-oriented deformation
        # field vs FS-native warp) before trusting this blindly.
        mri_warp_convert -g "${mni_t1}" \
            --inras "${xfm_dir}/mni_to_subj_warp.mgz" \
            --outm3z "${xfm_dir}/mni_to_subj.m3z"
    fi

    # --- Step 4: apply the warp to the Buckner atlas (nearest-neighbor, it's a label map) ---
    mri_vol2vol --mov "${buckner_atlas}" \
        --targ "${SUBJECTS_DIR}/${subj}/mri/norm.mgz" \
        --m3z "${xfm_dir}/mni_to_subj.m3z" \
        --noDefM3zPath \
        --interp nearest \
        --o "${buckner_fs}"

    # --- Step 5: restrict to cerebellum (aseg 8 = L-Cerebellum-Cortex, 47 = R-Cerebellum-Cortex) ---
    mri_binarize --i "${SUBJECTS_DIR}/${subj}/mri/aseg.mgz" \
        --match 8 47 \
        --o "${SUBJECTS_DIR}/${subj}/mri/cerebellum_mask.mgz"

    mri_mask "${buckner_fs}" \
        "${SUBJECTS_DIR}/${subj}/mri/cerebellum_mask.mgz" \
        "${buckner_fs}"

    rm -f "${SUBJECTS_DIR}/${subj}/mri/cerebellum_mask.mgz" \
          "${xfm_dir}/mni_in_subj_QC.mgz"
}


# Generate the lookup tables for a combined "BASE+Buckner" atlas:
#   atlas-${atlas}_Buckner_modified.txt : Buckner labels with IDs offset
#       past the base atlas' max node ID (input for labelconvert)
#   atlas-${atlas}_modified.txt        : base + offset Buckner labels
#       (input for header2timeseries.py, i.e. the header of the merged
#       timeseries output)
# The offset is computed dynamically instead of hardcoded, so it stays
# correct if the parcellation's max node count ever changes.
make_buckner_label_files() {
    local atlas=$1
    local base=${atlas%%+*}
    local anatdir="${fmriprepdir}/${subj}${sessionpath}anat"
    local pref="${subj}${sessionfile}"

    if [[ ! -f "${buckner_lut}" ]]; then
        echo -e "${RED}ERROR: Buckner lookup table not found at ${buckner_lut}; exiting${NC}" >&2
        exit 1
    fi

    local base_orig base_mod
    read -r base_orig base_mod < <(get_label_files "${base}")

    # Offset Buckner label IDs to avoid collision with the base atlas'
    # node IDs.
    local Nnodes
    Nnodes=$(awk 'BEGIN{m=0} $1+0>m{m=$1} END{print m}' "${base_mod}")
    local buckner_mod="${anatdir}/${pref}atlas-${atlas}_Buckner_modified.txt"
    awk -v c="${Nnodes}" -v OFS='\t' '{ $1 = $1 + c; print }' \
        "${buckner_lut}" > "${buckner_mod}"

    # combined table = base parcellation + offset Buckner labels
    cat "${base_mod}" "${buckner_mod}" \
        > "${anatdir}/${pref}atlas-${atlas}_modified.txt"
}

# Build the BOLD-space dseg of a combined "BASE+Buckner" atlas:
#   1. the base (Schaefer) parcellation via the standard bold-space path
#   2. the Buckner atlas resampled into BOLD space and relabelled with the
#      offset IDs (labelconvert)
#   3. both merged with fslmaths -add; after labelconvert the base dseg is
#      0 in the cerebellum and the Buckner dseg is 0 in the cortex, so a
#      plain add cannot create label collisions
build_buckner_combined_dseg() {
    local atlas=$1
    local base=${atlas%%+*}
    local anatdir="${fmriprepdir}/${subj}${sessionpath}anat"
    local pref="${subj}${sessionfile}"

    local base_dseg="${anatdir}/${pref}space-${outputspace}_atlas-${base}_dseg.nii.gz"
    local buckner_dseg="${anatdir}/${pref}space-${outputspace}_atlas-${atlas}_Buckner_dseg.nii.gz"
    local combined_dseg="${anatdir}/${pref}space-${outputspace}_atlas-${atlas}_dseg.nii.gz"

    echo -e "${BLUE}atlas = ${atlas} (Schaefer ${base} + Buckner cerebellum)${NC}"

    # subject-space Buckner atlas (no-op if already warped)
    warp_buckner_atlas

    # 1. cortical parcellation in BOLD space (standard Schaefer path)
    warp_atlas_to_bold_space "${base}"

    # 2. Buckner parcellation in BOLD space with offset label IDs
    make_buckner_label_files "${atlas}"
    local buckner_mod="${anatdir}/${pref}atlas-${atlas}_Buckner_modified.txt"


    echo -e "${BLUE}transform Buckner parcellation to boldspace${NC}"
    # follow the script's convention of a cached NIfTI copy (same as the
    # other +aseg atlases) before regridding to BOLD space
    local buckner_fs="${SUBJECTS_DIR}/${subj}/mri/Buckner2011_17Networks.mgz"
    local buckner_nii="${SUBJECTS_DIR}/${subj}/mri/Buckner2011_17Networks.nii.gz"
    [[ -f "${buckner_nii}" ]] || mrconvert "${buckner_fs}" "${buckner_nii}" -force
    local tempfile="${anatdir}/${pref}space-${outputspace}_atlas-${atlas}_Buckner_temp.nii.gz"
       
    mrgrid "${buckner_nii}" regrid - -template "${boldreffile}" -interp nearest \
    | mrconvert - "${tempfile}" -strides "$(mrinfo "${boldreffile}" -strides | tr -d '[]' | xargs | tr ' ' ',')"

    labelconvert "${tempfile}" "${buckner_lut}" "${buckner_mod}" "${buckner_dseg}" -force
    rm -f "${tempfile}"


    # 3. merge Buckner into the base parcellation
    echo -e "${BLUE}merge ${base} + Buckner${NC}"
    fslmaths "${base_dseg}" -add "${buckner_dseg}" "${combined_dseg}"

    rm -f "${buckner_dseg}" "${buckner_mod}"

}

# Resample a FreeSurfer-space atlas into the subject's BOLD space and
# relabel/sort it via labelconvert. Works for any atlas in ${ATLASES[@]}
# (Schaefer, BNA or combined "BASE+Buckner") since they all produce a
# subject-space parcellation volume.
warp_atlas_to_bold_space() {
    local atlas=$1
    local outfile="${fmriprepdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-${outputspace}_atlas-${atlas}_dseg.nii.gz"

    [[ -f "${outfile}" ]] && return 0

    # combined cortical + cerebellar ("BASE+Buckner") atlases are built from
    # the base parcellation dseg plus an offset Buckner dseg
    if [[ "${atlas}" == *+Buckner ]]; then
        build_buckner_combined_dseg "${atlas}"
        return 0
    fi

    echo -e "${YELLOW}atlas = ${atlas}${NC}"
    local parcfile="${SUBJECTS_DIR}/${subj}/mri/${atlas}+aseg.nii.gz"
    local tempfile="${fmriprepdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-${outputspace}_atlas-${atlas}_temp.nii.gz"

    echo -e "${BLUE}transform ${atlas} parcellation to boldspace${NC}"
    
  
    mrgrid "${parcfile}" regrid - -template "${boldreffile}" -interp nearest \
    | mrconvert - "${tempfile}" -strides "$(mrinfo "${boldreffile}" -strides | tr -d '[]' | xargs | tr ' ' ',')"


    local origfile modfile
    read -r origfile modfile < <(get_label_files "${atlas}")
    labelconvert "${tempfile}" "${origfile}" "${modfile}" "${outfile}" -force

    rm -f "${tempfile}"
}

compute_minmask() {
    base=$(remove_ext "${denoisedimg}")
    fslmaths "${denoisedimg}" -Tmin -thrP 25 -bin "${base}_minmask.nii.gz"
    fslmaths "${denoisedimg}" -mul "${base}_minmask.nii.gz" "${base}_minmasked.nii.gz"
}

# Extract ROI timeseries + ROI volumes (masked and min-masked) for one atlas
extract_atlas_timeseries() {
    local atlas=$1
    local ts_tsv="${fmriprepdir}/${subj}${sessionpath}timeseries/${subj}${sessionfile}task-rest_space-${outputspace}_atlas-${atlas}_rec-${denoise_protocol}_desc-timeseries_bold.tsv"
    local rv_txt="${fmriprepdir}/${subj}${sessionfile}task-rest_space-${outputspace}_atlas-${atlas}_desc-roivols_bold.txt"

    if [[ -f "${ts_tsv}" && -f "${rv_txt}" ]]; then
        echo -e "${GREEN}timeseries already extracted for ${atlas}${NC}"
        return 0
    fi

    local dsegfile="${fmriprepdir}/${subj}${sessionpath}anat/${subj}${sessionfile}space-${outputspace}_atlas-${atlas}_dseg.nii.gz"
    local tsdir="${fmriprepdir}/${subj}${sessionpath}timeseries"
    local rvdir="${fmriprepdir}/${subj}${sessionpath}roivols"

    echo -e "${BLUE}extract timeseries for ${atlas}${NC}"
    fslmeants -i "${denoisedimg}" --label="${dsegfile}" \
        > "${tempdir}/${subj}${sessionfile}atlas-${atlas}_timeseries.txt"
    fslmeants -i "${base}_minmasked.nii.gz" --label="${dsegfile}" \
        > "${tempdir}/${subj}${sessionfile}atlas-${atlas}_timeseries_minmasked.txt"

    echo -e "${BLUE}extract roi volumes for ${atlas}${NC}"
    fslstats -K "${dsegfile}" "${boldreffile}" -V | awk '{ print $2 }' \
        > "${rv_txt}"
    # NOTE: the original script accidentally reused the unmasked denoised
    # image here, so "_minmasked" roi volumes were identical to the
    # unmasked ones. Fixed to use the min-masked image.
    # fslstats -K "${dsegfile}" "${base}_minmasked.nii.gz" -V | awk '{ print $2 }' \
    #     > "${tempdir}/${subj}${sessionfile}atlas-${atlas}_roivols_minmasked.txt"

    rsync -a ${tempdir}/${subj}${sessionfile}atlas-${atlas}*.txt "${tsdir}/"

    local origfile modfile
    read -r origfile modfile < <(get_label_files "${atlas}")

    # add headers to the timeseries files, then drop the raw .txt versions
    for suffix in "" "_minmasked"; do
        local raw="${tsdir}/${subj}${sessionfile}atlas-${atlas}_timeseries${suffix}.txt"
        local times_tsv="${tsdir}/${subj}${sessionfile}task-rest_space-${outputspace}_atlas-${atlas}_rec-${denoise_protocol}_desc-timeseries${suffix}_bold.tsv"

        "${atlasdir}/header2timeseries.py" --timeseriesfile "${raw}" --atlasids "${modfile}" --outfile "${times_tsv}"
        if [[ -f "${times_tsv}" ]]; then rm -f "${raw}"; fi

    done
}

########################################
# MAIN
########################################

load_modules
select_subject
check_denoised_image

mkdir -p "${fmriprepdir}/${subj}${sessionpath}anat"
boldreffile="${fmriprepdir}/${subj}${sessionpath}func/${subj}${sessionfile}task-rest_space-${outputspace}_boldref.nii.gz"

ATLASES=("${SCHAEFER_ATLASES[@]}")
if [[ "${INCLUDE_BNA}" == true ]]; then
    warp_bna_atlas
    ATLASES+=(BNA)
fi

ensure_fsfiles
for atlas in "${SCHAEFER_ATLASES[@]}"; do
    warp_schaefer_atlas "${atlas}"
done

for atlas in "${ATLASES[@]}"; do
    warp_atlas_to_bold_space "${atlas}"
done

mkdir -p "${fmriprepdir}/${subj}${sessionpath}timeseries" "${fmriprepdir}/${subj}${sessionpath}roivols"
compute_minmask

for atlas in "${ATLASES[@]}"; do
    extract_atlas_timeseries "${atlas}"
done

echo "--------------------------"
echo -e "${GREEN}done with subject = ${subj}${NC}"
echo "--------------------------"
