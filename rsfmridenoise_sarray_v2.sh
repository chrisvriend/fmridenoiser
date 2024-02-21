#!/bin/bash

# (C) Chris Vriend - AmsUMC - 16-12-2022

#SBATCH --job-name=fmridenoise
#SBATCH --mem=16G
#SBATCH --partition=luna-cpu-short
#SBATCH --qos=anw-cpu
#SBATCH --cpus-per-task=1
#SBATCH --time=00-0:30:00
#SBATCH --nice=2000

# SLURM ARRAY SETTINGS - change the array length to maximally the number of subject folders in the
# freesurfer directory. the number behind the % represents the number of simultaneous processes.
# https://help.rc.ufl.edu/doc/SLURM_Job_Arrays
# If you want to change the number of simultaneous tasks of an active job, you can use scontrol:
# scontrol update ArrayTaskThrottle=<count> JobId=<jobID>

#SBATCH --array 1-118%6
#SBATCH -o fmridenoiser_%A_%a.log
#SBATCH --mail-type=END,FAIL

# By default in SLURM, the emails for events BEGIN, END and FAIL apply to the job array as a whole rather than individual tasks. So:
# If you want per task emails,
##- #SBATCH --mail-type=BEGIN,END,FAIL,ARRAY_TASKS

# usage instructions
Usage() {
    cat <<EOF


    (C) Chris Vriend - AmsUMC - 16-12-2022
    script to denoise resting-state functional MRI data processed using fmriprep
    using one of several denoising strategies (see PIPELINE OPTIONS)
    this denoised scan can be used for further processing (e.g. timeseries extraction, ICA or seed-based connectivity)
    a report folder is additionally created with a html file to show spatial maps of the denoising regressors
    Note some 'future' and 'Deprecation' warnings will be shown but this does not affect the results!

    Usage: sbatch ./rsfmridenoise_sbatch.sh <fmriprepdir> <session> 
    Obligatory:
    fmriprepdir = full path to fmriprep directory that contains subject's fmrirep output

    Optional:
	session = session ID of fmriprep output, e.g. ses-T0. keep empty if there are no sessions

	additional options and paths may need to be checked/modified in the script

EOF
    exit 1
}

[ _$3 = _ ] && Usage

######################
# to add some color to the output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m'
######################

headdir=${1}
subjectfile=${2}
denoise_protocol=${3} # denoising_protocol, e.g 6HMP8PhysSpikeRegGS
task=${4}
sess=${5} # select session
run=${6}  # run, run-X


subj=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${subjectfile})
# random delay
duration=$((RANDOM % 40 + 2))
echo -e "${YELLOW}INITIALIZING...(wait a sec)${NC}"
echo
sleep ${duration}

if [ -z "$sess" ]; then
    # sess empty
    sessionpath=/
    sessionfile=_
else
    sessionpath=/${sess}/
    sessionfile=_${sess}_
fi

if [ -z "$run" ]; then
    # sess empty
    runfile=_
else
    runfile=_${run}_
fi

scriptdir=/scratch/anw/share-np/fmridenoiser/denoiser-1.0.1
# load modules
module load fsl/6.0.6.5
module load Anaconda3/2023.03
synthstrip=/scratch/anw/share-np/fmridenoiser/synthstrip.1.2.sif
export APPTAINER_BINDPATH="/data/anw/anw-gold,/scratch/anw/"

# activate denoiser environment
conda activate /scratch/anw/share/python-env/denoiserenv

###################
# input variables
###################
BBTHRESH=10     # Brain background threshold
FWHM=6          # Smoothing kernel (FWSE)
Ndummy=3        # number of dummy scans to remove
lpfilter=0.08   # low pass filter | typical freq band =0.01 - 0.08 Hz
hpfilter=0.009  # highpass filter
outputspace=T1w # output space

#MNI152NLin6Asym_res-2
#T1w
#####################
# PIPELINE OPTIONS
#####################

#24HMP8PhysSpikeReg
#24HMP8PhysSpikeReg4GS
#24HMPaCompCorSpikeReg - not yet implementend
#24HMPaCompCorSpikeReg4GS - not yet implemented
#ICAAROMA8Phys
#ICAAROMA8Phys4GS
######################

############## - START PROCESSING - ##############
cd ${headdir}

if [ -d ${headdir}/${subj}${sessionpath}func ]; then

    # input filenames
    funcimage=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-preproc_bold.nii.gz
    mean_func=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_boldref.nii.gz

    # funcjson=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-preproc_bold.json
    # TR=$(jq -r '.RepetitionTime' ${fmriprepdir}/${subj}${sessionpath}func/${funcjson})
    regressorfile=${subj}${sessionfile}task-${task}${runfile}desc-confounds_timeseries.tsv

    # define output filenames
    funcimagenodummy=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_desc-preproc-dummy_bold.nii.gz
    funcimagesmooth=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_rec-susan${FWHM}mm_bold.nii.gz
    funcimagesmoothNR=${subj}${sessionfile}task-${task}${runfile}space-${outputspace}_rec-susan${FWHM}mm_desc-${denoise_protocol}_bold.nii.gz

    if [ ! -f ${headdir}/${subj}${sessionpath}func/${funcimagesmoothNR} ]; then

        echo
        echo -e "${GREEN}--------------${NC}"
        echo -e "${GREEN}subj = ${subj}${NC}"
        echo -e "${GREEN}--------------${NC}"
        echo

        if [[ ${headdir} == *"/data/anw/anw-archive"* ]]; then

            origdir=${headdir}
            headdir=~/my-scratch/denoise_work
            mkdir -p "${headdir}/${subj}${sessionpath}"
            (cd ${origdir}/${subj}${sessionpath} && find -type f -name "*task-${task}_space-${outputspace}*" -exec rsync -avR {} ${headdir}/${subj}${sessionpath} \;)
            (cd ${origdir}/${subj}${sessionpath} && find -type f -name "*task-${task}_desc-confounds*" -exec rsync -avR {} ${headdir}/${subj}${sessionpath} \;)

        fi
        cd ${headdir}/${subj}${sessionpath}func
        outputdir=${headdir}/${subj}${sessionpath}func

        ########################
        # remove dummy scans
        ########################

        Nvols=$(fslnvols ${funcimage})
        rvols=$(echo "${Nvols}-${Ndummy}" | bc -l)
        echo -e "${YELLOW}No. vols in original scan = ${Nvols}${NC}"
        lregfile=$(echo "${Ndummy}+2" | bc -l)

        if [ ! -f ${funcimagenodummy} ]; then
            echo -e "${YELLOW}...removing dummy scans...${NC}"
            fslroi ${funcimage} ${funcimagenodummy} ${Ndummy} -1
            Dvols=$(fslnvols ${funcimagenodummy})
            echo -e "${YELLOW}No. vols in modified scan = ${Dvols}${NC}"

        else
            echo -e "${YELLOW}...dummy scans already removed...${NC}"
            Dvols=$(fslnvols ${funcimagenodummy})
            echo -e "${YELLOW}No. vols in modified scan = ${Dvols}${NC}"

        fi

        # modify regressors file
        rm -f task-${task}${runfile}confounders-dummy.tsv
        cat ${regressorfile} | head -1 >task-${task}${runfile}confounders-dummy.tsv
        echo "$(tail -n +${lregfile} ${regressorfile})" >>task-${task}${runfile}confounders-dummy.tsv

        rows=$(cat task-${task}${runfile}confounders-dummy.tsv | wc -l)
        clines=$(echo "${rvols} + 1" | bc -l)

        if test ${rows} -ne ${clines} || test $(fslnvols ${funcimagenodummy}) -ne ${rvols}; then
            echo -e "${RED}ERROR! number of lines in regressor file does not match the number of functional volumes${NC}"
            exit
        fi

        unset rows clines rvols lregfile

        ######################################
        # PREPARE FOR SMOOTHING USING FSL SUSAN
        ######################################
        echo

        # Perform skullstrip on the mean_func data (=boldref)
        echo -e "${YELLOW}skullstrip mean functional image${NC}"
        echo
        apptainer run --cleanenv ${synthstrip} -i ${headdir}/${subj}${sessionpath}func/${mean_func} \
            -m ${headdir}/${subj}${sessionpath}func/mask.nii.gz
        echo
        # synthstrip does something weird to the header that leads to
        # warning messages in the next step. Therefore we clone the header
        # from the input image
        fslcpgeom ${headdir}/${subj}${sessionpath}func/${mean_func} \
            ${headdir}/${subj}${sessionpath}func/mask.nii.gz

        fslmaths ${headdir}/${subj}${sessionpath}func/${funcimagenodummy} \
            -mas ${headdir}/${subj}${sessionpath}func/mask.nii.gz \
            ${headdir}/${subj}${sessionpath}func/funcimage_BET

        echo -e "${YELLOW}calculate parameters for SUSAN...${NC}"
        lowerp=$(echo "scale=6; $(fslstats ${headdir}/${subj}${sessionpath}func/funcimage_BET -p 2)" | bc)
        upperp=$(echo "scale=6; $(fslstats ${headdir}/${subj}${sessionpath}func/funcimage_BET -p 98)" | bc)

        p98thr=$(echo "scale=6; ($upperp-$lowerp)/${BBTHRESH}" | bc)

        # Use fslmaths to threshold the brain extracted data based on the highpass filter above
        # use "mask" as a binary mask, and Tmin to specify we want the minimum across time
        fslmaths ${headdir}/${subj}${sessionpath}func/funcimage_BET -thr ${p98thr} \
            -Tmin -bin ${headdir}/${subj}${sessionpath}func/pre_thr_mask -odt char

        fslmaths ${headdir}/${subj}${sessionpath}func/funcimage_BET -mas ${headdir}/${subj}${sessionpath}func/pre_thr_mask \
            ${headdir}/${subj}${sessionpath}func/func_data_thresh
        # We now take this functional data , and create a "mean_func" image that is the mean across time (Tmean)
        fslmaths ${headdir}/${subj}${sessionpath}func/func_data_thresh -Tmean \
            ${headdir}/${subj}${sessionpath}func/mean_func

        # To run susan, FSLs tool for noise reduction, we need a brightness threshold.
        # Calculated using fslstats based on
        # https://neurostars.org/t/smoothing-images-by-susan-after-fmriprep/16453/4
        medint=$(fslstats ${headdir}/${subj}${sessionpath}func/funcimage_BET \
            -k ${headdir}/${subj}${sessionpath}func/pre_thr_mask -p 50)
        brightthresh=$(echo "scale=6; ((${medint}*0.75))" | bc)
        #echo "brightthreshold = ${brightthresh}"

        #####################
        # SUSAN SMOOTHING
        #####################

        ssize=$(echo "scale=11; ((${FWHM}/2.355))" | bc)

        if [ ! -f ${headdir}/${subj}${sessionpath}func/${funcimagesmooth} ]; then
            echo -e "${YELLOW}...running SUSAN${NC}"
            # susan uses nonlinear filtering to reduce noise
            # by only averaging a voxel with local voxels which have similar intensity
            susan ${headdir}/${subj}${sessionpath}func/func_data_thresh ${brightthresh} ${ssize} 3 1 1 \
                ${headdir}/${subj}${sessionpath}func/mean_func ${medint} \
                ${headdir}/${subj}${sessionpath}func/func_data_smooth
            # 3 means 3D smoothing
            # 1 says to use a local median filter
            # 1 says that we determine the smoothing area from 1 secondary image, "mean_func" and then we use the same brightness threshold for the secondary image.
            # prefiltered_func_data_smooth is the output image

            # Now we mask the smoothed functional data with the mask image, and overwrite the smoothed image.
            fslmaths ${headdir}/${subj}${sessionpath}func/func_data_smooth \
                -mas ${headdir}/${subj}${sessionpath}func/mask ${headdir}/${subj}${sessionpath}func/func_data_smooth

            mv ${headdir}/${subj}${sessionpath}func/func_data_smooth.nii.gz \
                ${headdir}/${subj}${sessionpath}func/${funcimagesmooth}

        fi

        cd ${headdir}/${subj}${sessionpath}func

        ###############################
        # DEFINE CONFOUND REGRESSORS
        ##############################
        echo
        echo -e "${BLUE}define confound regressors${NC}"
        echo
        if [[ ${denoise_protocol} == "ICAAROMA8Phys" ]] ||
            [[ ${denoise_protocol} == "ICAAROMA8Phys4GS" ]]; then

            aromaregressors=${subj}${sessionfile}task-${task}_AROMAnoiseICs.csv

            rm -f aroma_motion_regressors.txt
            OLDIFS=$IFS
            IFS=","

            for c in $(cat ${aromaregressors}); do

                if test $c -lt 10; then
                    cc=0${c}
                elif test $c -lt 100; then
                    cc=${c}
                fi

                echo aroma_motion_${cc} >>aroma_motion_regressors.txt

            done
            IFS=$OLDIFS

            # delim=""
            # for item in $(cat aroma_motion_regressors.txt); do
            #   printf "%s" "$delim$item" >>aroma_comma.txt
            #   delim=","
            # done

            AROMA=$(cat aroma_motion_regressors.txt)
            #rm aroma_comma.txt
            # determine despiking regressors

        elif [[ ${denoise_protocol} == *"SpikeReg"* ]]; then

            ${scriptdir}/find_motionoutliers.py -dir ${headdir} -subjid ${subj} -session ${sess} -task ${task}
            # -run

            # delim=""
            # for item in $(cat ${subj}${sessionfile}motion_outliers.txt); do
            #     printf "%s" "$delim$item" >>motion_comma.txt
            #     delim=","
            # done

            # MOTIONOUT=$(cat motion_comma.txt)
            # rm -f motion_comma.txt

            echo -e "${YELLOW}number of motion outliers for despiking  = $(cat ${subj}${sessionfile}task-${task}${runfile}motion_outliers.txt | wc -l)${NC}"

            MOTIONOUT=($(cat ${subj}${sessionfile}task-${task}${runfile}motion_outliers.txt))

        elif [[ ${denoise_protocol} == "24HMP8Phys" || ${denoise_protocol} == "24HMP8PhysGS" || ${denoise_protocol} == "24HMP8Phys4GS" || ${denoise_protocol} == "6HMP8PhysGS" || ${denoise_protocol} == "6HMP8Phys4GS" ]]; then
            :

        else
            echo -e "${RED}ERROR! denoise protocol not defined${NC}"
            exit
        fi

        if [[ ${denoise_protocol} == "ICAAROMA8Phys" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 ${AROMA}")
        elif [[ ${denoise_protocol} == "ICAAROMA8Phys4GS" ]]; then
            confounds=$(echo "global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2 csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 ${AROMA}")
        elif [[ ${denoise_protocol} == "24HMP8PhysSpikeReg" ]]; then
            if [ -z ${MOTIONOUT} ]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2")
            else
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 ${MOTIONOUT[@]}")
            fi
        elif [[ ${denoise_protocol} == "24HMP8PhysSpikeReg4GS" ]]; then
            if [ -z ${MOTIONOUT} ]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2")
            else
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2 ${MOTIONOUT[@]}")
            fi
        elif [[ ${denoise_protocol} == "24HMP8PhysSpikeRegGS" ]]; then
            if [ -z ${MOTIONOUT} ]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal")
            else
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal ${MOTIONOUT[@]}")
            fi
        elif [[ ${denoise_protocol} == "6HMP8PhysSpikeRegGS" ]]; then
            if [ -z ${MOTIONOUT} ]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal")
            else
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal ${MOTIONOUT[@]}")
            fi
        elif [[ ${denoise_protocol} == "6HMP8PhysSpikeReg4GS" ]]; then
            if [ -z ${MOTIONOUT} ]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2")
            else
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2 ${MOTIONOUT[@]}")
            fi
        elif [[ ${denoise_protocol} == "24HMP8Phys" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2")
        elif [[ ${denoise_protocol} == "24HMP8PhysGS" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal")
        elif [[ ${denoise_protocol} == "24HMP8Phys4GS" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2")
        elif [[ ${denoise_protocol} == "6HMP8PhysGS" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal")
        elif [[ ${denoise_protocol} == "6HMP8Phys4GS" ]]; then
            confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_y trans_z rot_x rot_y rot_z global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2")

        else

            echo -e "${RED}ERROR! pipeline for denoising not recognized${NC}"
            exit
        fi

        ###########################################################################
        # PERFORM DENOISING
        ###########################################################################

        if [ ! -f ${outputdir}/${funcimagesmoothNR} ]; then
            cd ${outputdir}
            echo
            echo -e "${YELLOW}denoising functional image${NC}"
            echo
            echo -e "${YELLOW}standby for some warning messages :)${NC}"
            sleep 2
            #cp ${scriptdir}/report_template.html ${headdir}/${subj}${sessionpath}func

            python ${scriptdir}/run_denoise.py ${outputdir}/${funcimagesmooth} \
                ${outputdir}/task-${task}${runfile}confounders-dummy.tsv \
                ${outputdir} \
                --lp_filter ${lpfilter} --hp_filter ${hpfilter} \
                --col_names ${confounds}
            sleep 1

            mv *bold_NR.nii.gz ${outputdir}/${funcimagesmoothNR}
            rm *bold_temp.nii.gz
            unset temp
        fi

        unset ssize uppert brightthresh lowerp upperp p98thr medmcf

        # print confounds to json sidecar
        confounds2=$(printf '%s\n' ${confounds})
        confoundarray=$(printf '%s\n' ${confounds2[@]} | jq -R . | jq -s .)
        temp=$(remove_ext ${funcimagesmoothNR})
        echo "
      {
        \"Confounds\": ${confoundarray},
        \"Smoothed\":${FWHM},
        \"IntendedFor\": \"${funcimagesmoothNR}\"
      }
      " >${outputdir}/${temp}.json
        temp2=$(remove_ext ${funcimage})
        #combine json files into one
        echo "$(jq -s 'add' ${outputdir}/${temp}.json ${headdir}/${subj}${sessionpath}func/${temp2}.json)" >${temp}.json
        #update taskname tag
        taskvar="${task}"
        echo "$(jq --arg urlarg "${taskvar}" '. += {"TaskName": $urlarg }' ${temp}.json)" >${temp}.json
        # update skull stripped tage
        echo "$(jq '.SkullStripped = true' ${temp}.json)" >${temp}.json
        unset temp
        # clean-up
        rm -f ${headdir}/${subj}${sessionpath}func/func_data_smooth_usan_size.nii.gz \
            ${headdir}/${subj}${sessionpath}func/funcimage_BET.nii.gz \
            ${headdir}/${subj}${sessionpath}func/mean_func.nii.gz \
            ${headdir}/${subj}${sessionpath}func/func_data_thresh.nii.gz \
            ${headdir}/${subj}${sessionpath}func/mask.nii.gz \
            ${headdir}/${subj}${sessionpath}func/pre_thr_mask.nii.gz
        # optional
        # rm ${headdir}/${subj}${sessionpath}func/${funcimagesmooth}
        if [[ ! -z ${origdir} ]]; then

            rsync -av --ignore-existing ${headdir}/${subj}${sessionpath}func/ ${origdir}/${subj}${sessionpath}func/
            if [[ -f ${origdir}/${subj}${sessionpath}func/${funcimagesmoothNR} ]]; then
                rm -r ~/my-stratch/denoise_work/${subj}${sessionpath}func/
                unset origdir
            fi

        fi

        echo -e "${GREEN}--------------------------${NC}"
        echo -e "${GREEN}denoising done for ${subj} - ${sess}${NC}"
        echo -e "${GREEN}--------------------------${NC}"

    else
        echo -e "${GREEN}denoising already performed for ${subj}- ${sess}${NC}"
        echo -e "${GREEN}--------------------------${NC}"

    fi # smdesnoised exist
else

    echo -e "${BLUE}--------------------------${NC}"
    echo -e "${BLUE}no fmriprep output found for ${subj} - ${sess}${NC}"
    echo -e "${BLUE}--------------------------${NC}"

fi
