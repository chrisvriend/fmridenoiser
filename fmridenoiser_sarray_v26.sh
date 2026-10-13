#!/bin/bash
#
# (C) Chris Vriend - AmsUMC - 16-12-2022
# rsfmridenoise - optimized/refactored version
# 12-10-2026: added SMOOTHTYPE option to choose between SUSAN smoothing
#             (default, original behaviour) and regular Gaussian FWHM
#             smoothing (fslmaths -kernel gauss -fmean)
# 12-10-2026: SMOOTHTYPE and outputspace are now user inputs, settable via
#             the --smoothtype and --space command-line options
#
#SBATCH --job-name=fmridenoise
#SBATCH --mem=4G
#SBATCH --partition=defq
#SBATCH --cpus-per-task=1
#SBATCH --time=00-0:30:00
#SBATCH --nice=2000
#SBATCH --array=1-603%40
#SBATCH --output=%x_%A_%a.log
#SBATCH --mail-type=END,FAIL

Usage() {
    cat <<EOF

    (C) Chris Vriend - AmsUMC - 09-10-2026
    script to denoise resting-state functional MRI data processed using fmriprep
    using one of several denoising strategies (see PIPELINE OPTIONS in the script)
    the denoised scan can be used for further processing (e.g. timeseries
    extraction, ICA or seed-based connectivity)
    a report folder is additionally created with a html file to show spatial maps
    of the denoising regressors
    Note some 'future' and 'Deprecation' warnings will be shown but this does not
    affect the results!

    Usage: sbatch ./fmridenoiser_sarray_v26.sh <fmriprepdir> [subjectfile] <denoise_protocol> <task> [session] [run] [options]

    Obligatory:
      fmriprepdir      = full path to fmriprep directory that contains subjects' fmriprep output
      denoise_protocol = denoising protocol, e.g. 6HMP8PhysSpikeRegGS (see PIPELINE OPTIONS)
      task             = task label, e.g. rest

    Optional:
      subjectfile      = text file with one subject ID per line (one line per array task)
                         if omitted, the subjects are collected automatically from the
                         sub-* directories in the fmriprep directory; make sure the
                         --array SBATCH option then covers the number of subjects
                         (array tasks beyond the number of subjects exit gracefully)
      session = session ID of fmriprep output, e.g. ses-T0 (keep empty if there are no sessions)
      run     = run label, e.g. run-1

    Options (user inputs; may be placed anywhere after the script name):
      --smoothtype <susan|gauss>
                         smoothing method (default: susan)
                           susan = FSL SUSAN nonlinear, edge-preserving smoothing
                           gauss = regular Gaussian FWHM smoothing (fslmaths -kernel gauss -fmean)
      --space <space>    output space of the fmriprep data (default: MNI152NLin6Asym_res-2)
                         e.g. T1w or MNI152NLin6Asym_res-2
      the smoothing method is recorded in the output filename (rec-susan<FWHM>mm vs
      rec-gauss<FWHM>mm) and in the json sidecar ("SmoothingMethod"); the output
      space is part of all input and output filenames (space-<space>)
      both options can also be set via the environment (SMOOTHTYPE / OUTPUTSPACE)
      or the CONFIG defaults; the command line takes precedence

    examples:
      sbatch ./fmridenoiser_sarray_v26.sh /path/to/fmriprep 6HMP8PhysSpikeRegGS rest
      sbatch ./fmridenoiser_sarray_v26.sh /path/to/fmriprep 6HMP8PhysSpikeRegGS rest --smoothtype gauss
      sbatch ./fmridenoiser_sarray_v26.sh /path/to/fmriprep subjects.txt 24HMP8PhysGS rest ses-T0 run-1 --space T1w --smoothtype gauss

    additional options and paths may need to be checked/modified in the CONFIG section

EOF
    exit 1
}

# ------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------
scriptdir=/appdata/users/P042819/git_fmridenoiser/denoiser-1.0.1

BBTHRESH=10      # Brain background threshold (SUSAN only)
FWHM=6           # Smoothing kernel (FWHM)
Ndummy=3         # number of dummy scans to remove
lpfilter=0.08    # low pass filter  | typical freq band = 0.01 - 0.08 Hz
hpfilter=0.009   # highpass filter

# user-settable options; precedence (highest first):
#   1. command line:  --smoothtype <susan|gauss>   and   --space <space>
#   2. environment:   SMOOTHTYPE=... sbatch ...   /   OUTPUTSPACE=... sbatch ...
#   3. the defaults below
# smoothing method:
#   susan = FSL SUSAN nonlinear, edge-preserving smoothing (original behaviour)
#   gauss = regular Gaussian FWHM smoothing (fslmaths -kernel gauss -fmean)
SMOOTHTYPE="${SMOOTHTYPE:-susan}"
# output space: T1w or MNI152NLin6Asym_res-2
outputspace="${OUTPUTSPACE:-MNI152NLin6Asym_res-2}"

# grand mean scaling of the final output (after denoising and smoothing)
# override per submission with e.g. DOGMS=false sbatch ...
DOGMS="${DOGMS:-true}"
GMSTARGET=10000

# PIPELINE OPTIONS (denoise_protocol):
#   6HMP8PhysGS | 6HMP8Phys4GS | 6HMP8PhysSpikeRegGS | 6HMP8PhysSpikeReg4GS
#   24HMP8Phys | 24HMP8PhysGS | 24HMP8Phys4GS
#   24HMP8PhysSpikeReg | 24HMP8PhysSpikeRegGS | 24HMP8PhysSpikeReg4GS
#   ICAAROMA8Phys | ICAAROMA8Phys4GS
#   24HMPaCompCorSpikeReg[4GS] - not yet implemented
# ------------------------------------------------------------------

# output colours
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'; BLUE='\033[0;36m'; NC='\033[0m'
err() { echo -e "${RED}ERROR: ${*}${NC}" >&2; exit 1; }

# ------------------------------------------------------------------
# command-line arguments
# ------------------------------------------------------------------
# optional named options are parsed first so that the positional arguments
# keep their original order and meaning:
#   --smoothtype <susan|gauss>   smoothing method (default: susan)
#   --space <space>              output space (default: MNI152NLin6Asym_res-2)
# the --option=value form is also accepted; unknown options are rejected so
# that typos fail immediately instead of silently changing the analysis
positional=()
while [ $# -gt 0 ]; do
    case "${1}" in
        --smoothtype)
            [ $# -ge 2 ] || err "--smoothtype requires a value: susan or gauss"
            SMOOTHTYPE="${2}"
            shift 2
            ;;
        --smoothtype=*)
            SMOOTHTYPE="${1#--smoothtype=}"
            shift
            ;;
        --space|--outputspace)
            [ $# -ge 2 ] || err "--space requires a value, e.g. T1w or MNI152NLin6Asym_res-2"
            outputspace="${2}"
            shift 2
            ;;
        --space=*|--outputspace=*)
            outputspace="${1#*=}"
            shift
            ;;
        -h|--help)
            Usage
            ;;
        -*)
            err "unknown option ${1} (see Usage)"
            ;;
        *)
            positional+=("${1}")
            shift
            ;;
    esac
done
set -- "${positional[@]}"

# minimum call: <fmriprepdir> <denoise_protocol> <task>
if [ $# -lt 3 ]; then
    Usage
fi

fmriprepdir=$1
shift

# subjectfile is optional:
#   - an explicitly passed empty string ("") is ignored, so the old call style
#     <fmriprepdir> "" <denoise_protocol> <task> keeps working
#   - if the next argument is an existing file it is used as the subject file
#   - otherwise it is assumed to be the denoise protocol and the subjects are
#     collected automatically from the fmriprep directory
if [ $# -gt 0 ] && [ -z "${1}" ]; then
    shift
fi

subjectfile=""
if [ $# -gt 0 ] && [ -n "${1}" ]; then
    if [ -f "${1}" ]; then
        subjectfile="${1}"
        shift
    elif [[ "${1}" == *.txt || "${1}" == */* ]]; then
        # looks like an intended subject file that does not exist
        err "subject file ${1} not found"
    fi
fi

denoise_protocol=${1:-}
task=${2:-}
sess=${3:-}
run=${4:-}

[ -n "${denoise_protocol}" ] || err "denoise_protocol not specified"
[ -n "${task}" ] || err "task not specified"

# validate the user inputs early so an invalid value fails fast
case "${SMOOTHTYPE}" in
    susan|gauss) ;;
    *) err "SMOOTHTYPE must be 'susan' or 'gauss' (got '${SMOOTHTYPE}')" ;;
esac

# extend this list if you denoise data in other fmriprep output spaces
case "${outputspace}" in
    T1w|MNI152NLin6Asym_res-2|MNI152NLin6Asym|MNI152NLin2009cAsym_res-2|MNI152NLin2009cAsym) ;;
    *) err "outputspace '${outputspace}' not recognized; expected e.g. T1w or MNI152NLin6Asym_res-2 (extend the list in the script if you use other spaces)" ;;
esac

# ------------------------------------------------------------------
# environment (conda hooks are not 'set -u' safe, so activate first)
# ------------------------------------------------------------------
module load fsl anaconda3 freesurfer || err "module load failed"
export APPTAINER_BINDPATH="/scratch/users/P042819,/net/beegfs/users/P042819/"

conda activate fmridenoiser || err "failed to activate conda environment fmridenoiser"

set -euo pipefail
trap 'echo -e "${RED}FAILED at line ${LINENO}${NC}" >&2' ERR

# ------------------------------------------------------------------
# subject selection and path/filename setup
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# subject list
#   either read from the subject file or collected automatically from the
#   sub-* directories in the fmriprep directory
# ------------------------------------------------------------------
if [ -n "${subjectfile}" ]; then
    [ -f "${subjectfile}" ] || err "subject file ${subjectfile} not found"
    mapfile -t subjects < <(tr -d '\r' < "${subjectfile}" | grep -v '^[[:space:]]*$' | sort -u)
    echo -e "${YELLOW}read ${#subjects[@]} subject(s) from ${subjectfile}${NC}"
else
    [ -d "${fmriprepdir}" ] || err "fmriprep directory ${fmriprepdir} not found"
    mapfile -t subjects < <(find "${fmriprepdir}" -mindepth 1 -maxdepth 1 -type d -name 'sub-*' -printf '%f\n' | sort)
    echo -e "${YELLOW}collected ${#subjects[@]} subject(s) from ${fmriprepdir}${NC}"
fi

if [ ${#subjects[@]} -eq 0 ]; then
    err "no subjects found"
fi

# SLURM_ARRAY_TASK_ID is 1-based here; bash arrays are 0-based
idx=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${idx}" -ge "${#subjects[@]}" ]; then
    echo -e "${BLUE}array task ${SLURM_ARRAY_TASK_ID} exceeds the number of subjects (${#subjects[@]}) - nothing to do${NC}"
    exit 0
fi
subj="${subjects[${idx}]}"

if [ -n "${sess}" ]; then
    sessionpath="/${sess}"; sessionfile="_${sess}_"
else
    sessionpath=""; sessionfile="_"
fi
if [ -n "${run}" ]; then
    runfile="_${run}_"
else
    runfile="_"
fi

funcdir="${fmriprepdir}/${subj}${sessionpath}/func"

funcimage="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-preproc_bold.nii.gz"
mean_func="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_boldref.nii.gz"
regressorfile="${subj}${sessionfile}task-${task}${runfile}desc-confounds_timeseries.tsv"
funcimagenodummy="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-preproc-dummy_bold.nii.gz"
funcimagemask="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-brain_mask.nii.gz"
# the smoothing method is encoded in the filenames so that susan and gauss
# outputs never overwrite each other and the resume-checks are per method
funcimagesmooth="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_rec-${SMOOTHTYPE}${FWHM}mm_bold.nii.gz"
funcimagesmoothNR="${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_rec-${SMOOTHTYPE}${FWHM}mm_desc-${denoise_protocol}_bold.nii.gz"
confoundersdummy="task-${task}${runfile}confounders-dummy.tsv"

# ------------------------------------------------------------------
# early exits
# ------------------------------------------------------------------
if [ ! -d "${funcdir}" ]; then
    echo -e "${BLUE}no fmriprep output found for ${subj} - ${sess}${NC}"
    exit 0
fi
if [ -f "${funcdir}/${funcimagesmoothNR}" ]; then
    echo -e "${GREEN}denoising already performed for ${subj} - ${sess}${NC}"
    exit 0
fi

# stagger array task start times to spread the I/O load
echo -e "${YELLOW}INITIALIZING...(wait a sec)${NC}"
sleep $((RANDOM % 40 + 2))

for f in "${funcimage}" "${mean_func}" "${regressorfile}"; do
    [ -f "${funcdir}/${f}" ] || err "missing input file ${funcdir}/${f}"
done

echo -e "${GREEN}--------------${NC}"
echo -e "${GREEN}subj = ${subj}${NC}"
echo -e "${GREEN}smoothing = ${SMOOTHTYPE} (FWHM = ${FWHM} mm)${NC}"
echo -e "${GREEN}outputspace = ${outputspace}${NC}"
echo -e "${GREEN}--------------${NC}"


# ------------------------------------------------------------------
# if the data live on the archive, copy the required files to scratch
# ------------------------------------------------------------------
origdir=""
if [[ "${fmriprepdir}" == *"/data/anw/anw-archive"* ]]; then
    origdir="${fmriprepdir}"
    fmriprepdir="${HOME}/my-scratch/denoise_work"
    funcdir="${fmriprepdir}/${subj}${sessionpath}/func"
    mkdir -p "${funcdir}"
    ( cd "${origdir}/${subj}${sessionpath}" &&
      find . -type f \( -name "*task-${task}_space-${outputspace}*" \
                       -o -name "*task-${task}_desc-confounds*" \) \
          -exec rsync -aR {} "${fmriprepdir}/${subj}${sessionpath}/" + )
fi

# ------------------------------------------------------------------
# define confound regressors (done early so an invalid protocol fails fast)
# ------------------------------------------------------------------
echo -e "${BLUE}define confound regressors${NC}"

# building blocks
PHYS8="csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2"
MOT6="trans_x trans_y trans_z rot_x rot_y rot_z"
MOT24="trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_power2 trans_y_derivative1_power2 trans_z trans_z_derivative1 trans_z_power2 trans_z_derivative1_power2 rot_x rot_x_derivative1 rot_x_power2 rot_x_derivative1_power2 rot_y rot_y_derivative1 rot_y_power2 rot_y_derivative1_power2 rot_z rot_z_derivative1 rot_z_power2 rot_z_derivative1_power2"
GS1="global_signal"
GS4="global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2"

AROMA=""
MOTIONOUT=()

case "${denoise_protocol}" in
    ICAAROMA*)
        # AROMA noise components -> regressor names (zero-padded to 2 digits)
        aromaregressors="${funcdir}/${subj}${sessionfile}task-${task}_AROMAnoiseICs.csv"
        [ -f "${aromaregressors}" ] || err "missing ${aromaregressors}"
        AROMA=$(tr ',' '\n' < "${aromaregressors}" | awk 'NF { printf "aroma_motion_%02d ", $1 }')
        ;;
    *SpikeReg*)
        if [ -n "${sess}" ]; then
        # shellcheck disable=SC2086  # unquoted on purpose: preserves original behaviour when no session is given
          ${scriptdir}/find_motionoutliers.py -dir "${fmriprepdir}" -subjid "${subj}" -session ${sess} -task "${task}"
        # -run
        else
         ${scriptdir}/find_motionoutliers.py -dir "${fmriprepdir}" -subjid "${subj}" -task "${task}"
        fi
        motionfile="${funcdir}/${subj}${sessionfile}task-${task}${runfile}motion_outliers.txt"
        [ -f "${motionfile}" ] || err "missing ${motionfile}"
        mapfile -t MOTIONOUT < "${motionfile}"
        echo -e "${YELLOW}number of motion outliers for despiking = ${#MOTIONOUT[@]}${NC}"
        ;;
esac

MOTSTR=""
if [ ${#MOTIONOUT[@]} -gt 0 ]; then
    MOTSTR="${MOTIONOUT[*]}"
fi

case "${denoise_protocol}" in
    ICAAROMA8Phys)         confounds="${PHYS8} ${AROMA}" ;;
    ICAAROMA8Phys4GS)      confounds="${PHYS8} ${GS4} ${AROMA}" ;;
    24HMP8Phys)            confounds="${PHYS8} ${MOT24}" ;;
    24HMP8PhysGS)          confounds="${PHYS8} ${MOT24} ${GS1}" ;;
    24HMP8Phys4GS)         confounds="${PHYS8} ${MOT24} ${GS4}" ;;
    24HMP8PhysSpikeReg)    confounds="${PHYS8} ${MOT24} ${MOTSTR}" ;;
    24HMP8PhysSpikeRegGS)  confounds="${PHYS8} ${MOT24} ${GS1} ${MOTSTR}" ;;
    24HMP8PhysSpikeReg4GS) confounds="${PHYS8} ${MOT24} ${GS4} ${MOTSTR}" ;;
    6HMP8PhysGS)           confounds="${PHYS8} ${MOT6} ${GS1}" ;;
    6HMP8Phys4GS)          confounds="${PHYS8} ${MOT6} ${GS4}" ;;
    6HMP8PhysSpikeRegGS)   confounds="${PHYS8} ${MOT6} ${GS1} ${MOTSTR}" ;;
    6HMP8PhysSpikeReg4GS)  confounds="${PHYS8} ${MOT6} ${GS4} ${MOTSTR}" ;;
    *) err "denoise protocol ${denoise_protocol} not recognized (see PIPELINE OPTIONS; if you meant this to be a subject file, check that it exists)" ;;
esac

# ------------------------------------------------------------------
# remove dummy scans and align the confounds file with them
# ------------------------------------------------------------------
Nvols=$(fslnvols "${funcdir}/${funcimage}")
rvols=$((Nvols - Ndummy))
lregfile=$((Ndummy + 2))
echo -e "${YELLOW}No. vols in original scan = ${Nvols}${NC}"

if [ ! -f "${funcdir}/${funcimagenodummy}" ]; then
    echo -e "${YELLOW}...removing dummy scans...${NC}"
    fslroi "${funcdir}/${funcimage}" "${funcdir}/${funcimagenodummy}" "${Ndummy}" -1
else
    echo -e "${YELLOW}...dummy scans already removed...${NC}"
fi
Dvols=$(fslnvols "${funcdir}/${funcimagenodummy}")
echo -e "${YELLOW}No. vols in modified scan = ${Dvols}${NC}"

# keep the header, drop the first Ndummy regressor rows
awk -v start="${lregfile}" 'NR==1 || NR>=start' \
    "${funcdir}/${regressorfile}" > "${funcdir}/${confoundersdummy}"
rows=$(wc -l < "${funcdir}/${confoundersdummy}")

if [ "${rows}" -ne $((rvols + 1)) ] || [ "${Dvols}" -ne "${rvols}" ]; then
    err "number of lines in regressor file does not match the number of functional volumes"
fi

# ------------------------------------------------------------------
# skullstrip and smoothing
# ------------------------------------------------------------------
echo -e "${YELLOW}skullstrip mean functional image${NC}"
#mri_synthstrip -i "${funcdir}/${mean_func}" -m "${funcdir}/mask.nii.gz"
# synthstrip does something weird to the header that leads to warning messages in
# the next step; therefore clone the header from the input image
#fslcpgeom "${funcdir}/${mean_func}" "${funcdir}/mask.nii.gz"
python3 - "${funcdir}/${funcimagemask}" "$funcdir" <<'EOF'
import sys, os
import nibabel as nib, numpy as np

infile, outdir = sys.argv[1], sys.argv[2]
outfile = os.path.join(outdir, 'mask_clean.nii.gz')

if os.path.abspath(infile) == os.path.abspath(outfile):
    sys.exit("Refusing to overwrite the input file; choose a different output directory.")

img = nib.load(infile)
d = np.asanyarray(img.dataobj).astype(np.uint8)
out = nib.Nifti1Image(d, img.affine)
out.header.set_data_dtype(np.uint8)
out.header.set_slope_inter(1, 0)
nib.save(out, outfile)
print("Saved:", outfile)
EOF

fslmaths "${funcdir}/${funcimagenodummy}" -mas "${funcdir}/mask_clean.nii.gz" "${funcdir}/funcimage_BET"

if [ "${SMOOTHTYPE}" = "susan" ]; then
    # --------------------------------------------------------------
    # SUSAN smoothing (nonlinear, edge-preserving)
    # --------------------------------------------------------------
    echo -e "${YELLOW}calculate parameters for SUSAN...${NC}"
    read -r lowerp upperp < <(fslstats "${funcdir}/funcimage_BET" -p 2 -p 98)
    p98thr=$(echo "scale=6; (${upperp} - ${lowerp}) / ${BBTHRESH}" | bc)

    fslmaths "${funcdir}/funcimage_BET" -thr "${p98thr}" -Tmin -bin "${funcdir}/pre_thr_mask" -odt char
    fslmaths "${funcdir}/funcimage_BET" -mas "${funcdir}/pre_thr_mask" "${funcdir}/func_data_thresh"
    fslmaths "${funcdir}/func_data_thresh" -Tmean "${funcdir}/mean_func"

    medint=$(fslstats "${funcdir}/funcimage_BET" -k "${funcdir}/pre_thr_mask" -p 50)
    brightthresh=$(echo "scale=6; ${medint} * 0.75" | bc)

    ssize=$(echo "scale=11; ${FWHM} / 2.355" | bc)

    if [ ! -f "${funcdir}/${funcimagesmooth}" ]; then
        echo -e "${YELLOW}...running SUSAN${NC}"
        # susan uses nonlinear filtering to reduce noise by only averaging a voxel with
        # local voxels of similar intensity
        # (3 = 3D smoothing, 1 = use a local median filter,
        #  1 = derive the smoothing area from 1 secondary image "mean_func" using the
        #      same brightness threshold)
        susan "${funcdir}/func_data_thresh" "${brightthresh}" "${ssize}" 3 1 1 \
              "${funcdir}/mean_func" "${medint}" "${funcdir}/func_data_smooth"
        fslmaths "${funcdir}/func_data_smooth" -mas "${funcdir}/mask_clean.nii.gz" "${funcdir}/func_data_smooth"
        mv -f "${funcdir}/func_data_smooth.nii.gz" "${funcdir}/${funcimagesmooth}"
    fi
else
    # --------------------------------------------------------------
    # regular Gaussian FWHM smoothing (linear)
    # sigma = FWHM / (2 * sqrt(2 * ln(2))) = FWHM / 2.3548
    # note: the input is already masked, so voxels at the brain edge are
    # averaged with zeros outside the mask and become slightly darker; this is
    # inherent to Gaussian smoothing of masked data (the same happens in the
    # 'smooth after denoising' approach of the older script). The mask is
    # re-applied afterwards to zero out non-brain voxels.
    # --------------------------------------------------------------
    ssize=$(echo "scale=11; ${FWHM} / 2.3548" | bc)

    if [ ! -f "${funcdir}/${funcimagesmooth}" ]; then
        echo -e "${YELLOW}...running Gaussian smoothing (FWHM = ${FWHM} mm)${NC}"
        fslmaths "${funcdir}/funcimage_BET" -kernel gauss "${ssize}" -fmean "${funcdir}/func_data_smooth"
        fslmaths "${funcdir}/func_data_smooth" -mas "${funcdir}/mask_clean.nii.gz" "${funcdir}/func_data_smooth"
        mv -f "${funcdir}/func_data_smooth.nii.gz" "${funcdir}/${funcimagesmooth}"
    fi
fi

# ------------------------------------------------------------------
# denoising
# ------------------------------------------------------------------
if [ ! -f "${funcdir}/${funcimagesmoothNR}" ]; then
    cd "${funcdir}"
    echo -e "${YELLOW}denoising functional image${NC}"
    # shellcheck disable=SC2086  # intentional word splitting: one argument per regressor
    python "${scriptdir}/run_denoise.py" "${funcdir}/${funcimagesmooth}" \
        "${funcdir}/${confoundersdummy}" \
        "${funcdir}" \
        --lp_filter "${lpfilter}" --hp_filter "${hpfilter}" \
        --col_names ${confounds}

    # rename the denoised output to its final name
    shopt -s nullglob
    NRfiles=("${funcdir}"/*bold_NR.nii.gz)
    shopt -u nullglob
    [ ${#NRfiles[@]} -eq 1 ] || err "expected exactly one *bold_NR.nii.gz output, found ${#NRfiles[@]}"
    mv -f "${NRfiles[0]}" "${funcdir}/${funcimagesmoothNR}"
    rm -f "${funcdir}"/*bold_temp.nii.gz
fi

# ------------------------------------------------------------------
# grand mean scaling (optional)
#   applied to the FINAL image: every upstream step is linear or
#   scale-equivariant, so a global scalar commutes with all of them;
#   placing it last keeps QC in native units and guarantees the
#   delivered file has global mean = GMSTARGET
# ------------------------------------------------------------------
gmsmean=""
if [ "${DOGMS}" = "true" ]; then
    gmsmean=$(fslstats "${funcdir}/${funcimagesmoothNR}" -M | awk '{print $1}')
    [ -n "${gmsmean}" ] || err "could not compute global mean of ${funcimagesmoothNR}"
    # sanity check: must be a positive number (also guards against nan/0 -> divide by zero)
    awk -v m="${gmsmean}" 'BEGIN{exit !(m+0 > 0)}' \
        || err "invalid global mean '${gmsmean}' for ${funcimagesmoothNR}"

    fslmaths "${funcdir}/${funcimagesmoothNR}" -div "${gmsmean}" -mul "${GMSTARGET}" \
             "${funcdir}/${funcimagesmoothNR}" -odt float
    echo -e "${YELLOW}grand mean scaling: ${gmsmean} -> ${GMSTARGET}${NC}"
fi


# ------------------------------------------------------------------
# JSON sidecar
# ------------------------------------------------------------------
origbase="${funcimage%.nii.gz}"
[ -f "${funcdir}/${origbase}.json" ] || err "missing fmriprep sidecar ${funcdir}/${origbase}.json"

confoundarray=$(jq -n --arg c "${confounds}" '$c | split(" ") | map(select(length > 0))')
jq -n \
    --argjson conf "${confoundarray}" \
    --argjson smoothed "${FWHM}" \
    --arg smoothmethod "${SMOOTHTYPE}" \
    --arg gmsmean "${gmsmean}" \
    --argjson gmstarget "${GMSTARGET}" \
    --arg intended "${funcimagesmoothNR}" \
    --arg task "${task}" \
    --slurpfile orig "${funcdir}/${origbase}.json" \
    '$orig[0] + {
        Confounds: $conf,
        Smoothed: $smoothed,
        SmoothingMethod: $smoothmethod,
        GrandMeanScaling: (if $gmsmean == "" then false
                           else {Target: $gmstarget, OriginalGlobalMean: ($gmsmean | tonumber)} end),
        IntendedFor: $intended,
        TaskName: $task,
        SkullStripped: true
    }' > "${funcdir}/${funcimagesmoothNR%.nii.gz}.json"

# ------------------------------------------------------------------
# clean-up
# ------------------------------------------------------------------
rm -f "${funcdir}"/{func_data_smooth_usan_size,funcimage_BET,mean_func,func_data_thresh,mask_clean,pre_thr_mask}.nii.gz
# optional:
# rm -f "${funcdir}/${funcimagesmooth}"

# # copy results back to the archive and remove the scratch copy
# if [ -n "${origdir}" ]; then
#     rsync -a --ignore-existing "${funcdir}/" "${origdir}/${subj}${sessionpath}/func/"
#     if [ -f "${origdir}/${subj}${sessionpath}/func/${funcimagesmoothNR}" ]; then
#         rm -rf "${HOME}/my-scratch/denoise_work/${subj}${sessionpath}/func"
#     fi
# fi

echo -e "${GREEN}--------------------------${NC}"
echo -e "${GREEN}denoising done for ${subj} - ${sess}${NC}"
echo -e "${GREEN}--------------------------${NC}"
