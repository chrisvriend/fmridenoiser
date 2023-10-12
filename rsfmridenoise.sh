#!/bin/bash

# (C) Chris Vriend - AmsUMC - 16-12-2022
# script to denoise resting-state functional MRI data processed using fmriprep
# using one of several denoising strategies (see PIPELINE OPTIONS)
# the folder "denoised" is created in the func folder of the fmriprep output
# the denoised and smoothed output = ${subj}_${session}_task-rest_space-${outputspace}_desc-preproc_bold-smooth_NR.nii.gz
# this denoised scan can be used for further processing
# a report folder is additionally created with a html file to show spatial maps of the denoising variables

# note some 'future' and 'Deprecation' warnings will be shown but this does not affect the results!

headdir=/home/anw/cvriend/my-scratch/TIPICCO

scriptdir=/scratch/anw/share-np/fmridenoiser/denoiser-1.0.1
# load modules
module load fsl/6.0.6.5
module load Anaconda3/2023.03
synthstrip=/scratch/anw/share-np/fmridenoiser/synthstrip.1.2.sif

export APPTAINER_BINDPATH="/data/anw/anw-gold,/scratch/anw"

# activate denoiser environment
conda activate /scratch/anw/share/python-env/denoiserenv
###################
# input variables
###################
BBTHRESH=10                       # Brain background threshold
SKERN=6.0                         # Smoothing kernel (FWSE)
TR=2.2                            # repetition time of fMRI
Ndummy=3                          # number of dummy scans to remove
lpfilter=0.08                     # low pass filter | typical freq band =0.01 - 0.08 Hz
hpfilter=0.009                    # highpass filter
outputspace=MNI152NLin6Asym_res-2 # output space

######################
# to add some color to the output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;36m'
NC='\033[0m'
######################

#####################
# PIPELINE OPTIONS
#####################
denoise_protocol=ICAAROMA8Phys

#24HMP8PhysSpikeReg
#24HMP8PhysSpikeReg4GS
#24HMPaCompCorSpikeReg - not yet implementend
#24HMPaCompCorSpikeReg4GS - not yet implemented
#ICAAROMA8Phys
#ICAAROMA8Phys4GS
######################

# select session

session=ses-T0

############## - START PROCESSING - ##############
cd ${headdir}
for subj in $(ls -d sub-TIPICCO080); do

    if [ -d ${headdir}/${subj}/${session}/func ]; then

        # input filenames
        funcimage=${subj}_${session}_task-rest_space-${outputspace}_desc-preproc_bold.nii.gz
        mean_func=${subj}_${session}_task-rest_space-${outputspace}_boldref.nii.gz
        regressorfile=${subj}_${session}_task-rest_desc-confounds_timeseries.tsv

        # define output filenames
        funcimagenodummy=${subj}_${session}_task-rest_space-${outputspace}_desc-preproc_dummy_bold.nii.gz
        funcimagesmooth=${subj}_${session}_task-rest_space-${outputspace}_desc-smooth_bold.nii.gz
        funcimageNR=${subj}_${session}_task-rest_space-${outputspace}_desc-smooth_${denoise_protocol}_bold.nii.gz

        outputdir=${headdir}/${subj}/${session}/func/denoised

        if [ ! -f ${headdir}/${subj}/${session}/func/denoised/${funcimageNR} ]; then
            echo
            echo -e "${GREEN}--------------${NC}"
            echo -e "${GREEN}subj = ${subj}${NC}"
            echo -e "${GREEN}--------------${NC}"
            echo
            cd ${headdir}/${subj}/${session}/func

            mkdir -p ${outputdir}

            ########################
            # remove dummy scans
            ########################

            Nvols=$(fslnvols ${funcimage})
            rvols=$(echo "${Nvols}-${Ndummy}" | bc -l)
            echo -e "${YELLOW}No. vols in original scan = ${Nvols}${NC}"
            lregfile=$(echo "${Ndummy}+2" | bc -l)

            if [ ! -f ${funcimagenodummy} ]; then
                echo -e "${BLUE}...removing dummy scans...${NC}"
                fslroi ${funcimage} ${funcimagenodummy} ${Ndummy} -1
                Dvols=$(fslnvols ${funcimagenodummy})
                echo -e "${YELLOW}No. vols in modified scan = ${Dvols}${NC}"
            fi

            # modify regressors file
            rm -f confounders-dummy.tsv
            cat ${regressorfile} | head -1 >confounders-dummy.tsv
            echo "$(tail -n +${lregfile} ${regressorfile})" >>confounders-dummy.tsv

            rows=$(cat confounders-dummy.tsv | wc -l)
            clines=$(echo "${rvols} + 1" | bc -l)

            if test ${rows} -ne ${clines} || test $(fslnvols ${funcimagenodummy}) -ne ${rvols}; then
                echo -e "${RED}number of lines in regressor file does not match the number of functional volumes${NC}"
                exit
            fi

            unset rows clines rvols lregfile

            ######################################
            # PREPARE FOR SMOOTHING USING FSL SUSAN
            ######################################
            echo

            # Perform skullstrip on the mean_func data (=boldref)
            echo -e "${BLUE}skullstrip mean functional image${NC}"
            echo
            ${synthstrip} -i ${headdir}/${subj}/${session}/func/${mean_func} \
                -m ${headdir}/${subj}/${session}/func/denoised/mask.nii.gz
            echo
            # synthstrip does something weird to the header that leads to
            # warning messages in the next step. Therefore we clone the header
            # from the input image
            fslcpgeom ${headdir}/${subj}/${session}/func/${mean_func} \
                ${headdir}/${subj}/${session}/func/denoised/mask.nii.gz

            fslmaths ${headdir}/${subj}/${session}/func/${funcimagenodummy} \
                -mas ${headdir}/${subj}/${session}/func/denoised/mask.nii.gz \
                ${headdir}/${subj}/${session}/func/funcimage_BET

            echo -e "${BLUE}calculate parameters for SUSAN...${NC}"
            lowerp=$(echo "scale=6; $(fslstats ${headdir}/${subj}/${session}/func/funcimage_BET -p 2)" | bc)
            upperp=$(echo "scale=6; $(fslstats ${headdir}/${subj}/${session}/func/funcimage_BET -p 98)" | bc)

            p98thr=$(echo "scale=6; ($upperp-$lowerp)/${BBTHRESH}" | bc)

            # Use fslmaths to threshold the brain extracted data based on the highpass filter above
            # use "mask" as a binary mask, and Tmin to specify we want the minimum across time
            fslmaths ${headdir}/${subj}/${session}/func/funcimage_BET -thr ${p98thr} \
                -Tmin -bin ${headdir}/${subj}/${session}/func/denoised/pre_thr_mask -odt char

            fslmaths ${headdir}/${subj}/${session}/func/funcimage_BET -mas ${headdir}/${subj}/${session}/func/denoised/pre_thr_mask \
                ${headdir}/${subj}/${session}/func/denoised/func_data_thresh
            # We now take this functional data , and create a "mean_func" image that is the mean across time (Tmean)
            fslmaths ${headdir}/${subj}/${session}/func/denoised/func_data_thresh -Tmean \
                ${headdir}/${subj}/${session}/func/denoised/mean_func

            # To run susan, FSLs tool for noise reduction, we need a brightness threshold.
            # Calculated using fslstats based on
            # https://neurostars.org/t/smoothing-images-by-susan-after-fmriprep/16453/4
            medint=$(fslstats ${headdir}/${subj}/${session}/func/funcimage_BET \
                -k ${headdir}/${subj}/${session}/func/denoised/pre_thr_mask -p 50)
            brightthresh=$(echo "scale=6; ((${medint}*0.75))" | bc)
            #echo "brightthreshold = ${brightthresh}"

            #####################
            # SUSAN SMOOTHING
            #####################

            ssize=$(echo "scale=11; ((${SKERN}/2.355))" | bc)

            if [ ! -f ${headdir}/${subj}/${session}/func/denoised/${funcimagesmooth} ]; then
                echo -e "${BLUE}...running SUSAN${NC}"
                # susan uses nonlinear filtering to reduce noise
                # by only averaging a voxel with local voxels which have similar intensity
                susan ${headdir}/${subj}/${session}/func/denoised/func_data_thresh ${brightthresh} ${ssize} 3 1 1 \
                    ${headdir}/${subj}/${session}/func/denoised/mean_func ${medint} \
                    ${headdir}/${subj}/${session}/func/denoised/func_data_smooth
                # 3 means 3D smoothing
                # 1 says to use a local median filter
                # 1 says that we determine the smoothing area from 1 secondary image, "mean_func" and then we use the same brightness threshold for the secondary image.
                # prefiltered_func_data_smooth is the output image

                # Now we mask the smoothed functional data with the mask image, and overwrite the smoothed image.
                fslmaths ${headdir}/${subj}/${session}/func/denoised/func_data_smooth \
                    -mas ${headdir}/${subj}/${session}/func/denoised/mask ${headdir}/${subj}/${session}/func/denoised/func_data_smooth

                mv ${headdir}/${subj}/${session}/func/denoised/func_data_smooth.nii.gz \
                    ${headdir}/${subj}/${session}/func/denoised/${funcimagesmooth}

            fi

            cd ${headdir}/${subj}/${session}/func

            ###############################
            # DEFINE CONFOUND REGRESSORS
            ##############################
            echo
            echo -e "${BLUE}define confound regressors of protocol: ${denoise_protocol}${NC}"
            echo
            if [[ ${denoise_protocol} == "ICAAROMA8Phys" ]] ||
                [[ ${denoise_protocol} == "ICAAROMA8Phys4GS" ]]; then

                aromaregressors=${subj}_${session}_task-rest_AROMAnoiseICs.csv

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

            elif [[ ${denoise_protocol} == "24HMP8PhysSpikeReg" ]] ||
                [[ ${denoise_protocol} == "24HMP8PhysSpikeReg4GS" ]]; then

                ${scriptdir}/find_motionoutliers.py -dir ${fmriprepdir} -subjid ${subj} -session ${session}

                delim=""
                for item in $(cat ${subj}_${session}_motion_outliers.txt); do
                    printf "%s" "$delim$item" >>motion_comma.txt
                    delim=","
                done

                MOTIONOUT=$(cat motion_comma.txt)
                rm -f motion_comma.txt

                echo -e "${YELLOW}number of motion outliers for despiking  = $(cat ${subj}_${session}_motion_outliers.txt | wc -l)${NC}"

            else
                echo -e "${RED}denoise protocol not defined${NC}"
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
                    confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 ${MOTIONOUT}")
                fi
            elif [[ ${denoise_protocol} == "24HMP8PhysSpikeReg4GS" ]]; then
                if [ -z ${MOTIONOUT} ]; then
                    confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2")
                else
                    confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 global_signal global_signal_derivative1 global_signal_power2 global_signal_derivative1_power2 ${MOTIONOUT}")
                fi
            elif [[ ${denoise_protocol} == "24HMP8Phys" ]]; then
                confounds=$(echo "csf csf_derivative1 csf_derivative1_power2 csf_power2 white_matter white_matter_derivative1 white_matter_power2 white_matter_derivative1_power2 trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2")

                # WORK IN PROGRESS!!
                #  https://neurostars.org/t/fmriprep-acompcor-from-wm-csf-masks/19535/3
                #confounds=$(echo "trans_x trans_x_derivative1 trans_x_power2 trans_x_derivative1_power2 trans_y trans_y_derivative1 trans_y_derivative1_power2 trans_y_power2 trans_z trans_z_derivative1 trans_z_derivative1_power2 trans_z_power2 rot_x rot_x_derivative1 rot_x_derivative1_power2 rot_x_power2 rot_y rot_y_derivative1 rot_y_derivative1_power2 rot_y_power2 rot_z rot_z_derivative1 rot_z_derivative1_power2 rot_z_power2 ${MOTIONOUT}")
            else

                echo -e "${RED}pipeline for denoising not recognized${NC}"
                exit
            fi

            ###########################################################################
            # PERFORM DENOISING
            ###########################################################################

            if [ ! -f ${outputdir}/${funcimageNR} ]; then
                mkdir -p ${outputdir}/report
                cd ${headdir}/${subj}/${session}/func/denoised
                echo
                echo -e "${BLUE}denoising functional image${NC}"
                echo
                echo -e "${YELLOW}stand by for some warning messages :)${NC}"
                sleep 2
                cp ${scriptdir}/report_template.html ${headdir}/${subj}/${session}/func/denoised

                python ${scriptdir}/run_denoise.py ${funcimagesmooth} \
                    ${headdir}/${subj}/${session}/func/confounders-dummy.tsv \
                    ${outputdir} --out_figure_path ${outputdir}/report \
                    --lp_filter ${lpfilter} --hp_filter ${hpfilter} \
                    --col_names ${confounds}

                temp=$(remove_ext ${funcimagesmooth})
                mv ${temp}_NR.nii.gz ${funcimageNR}
                unset temp
            fi

            unset ssize uppert brightthresh lowerp upperp p98thr medmcf

            # print confounds to json sidecar
            confounds2=$(printf '%s\n' ${confounds})
            confoundarray=$(printf '%s\n' ${confounds2[@]} | jq -R . | jq -s .)
            temp=$(remove_ext ${funcimageNR})
            echo "
      {
        \"Confounds\": ${confoundarray},
        \"Smoothed\":true,
        \"IntendedFor\": \"${funcimageNR}\"
      }
      " >${outputdir}/${temp}.json
            temp2=$(remove_ext ${funcimage})
            #combine json files into one
            echo "$(jq -s 'add' ${outputdir}/${temp}.json ${headdir}/${subj}/${session}/func/${temp2}.json)" >${temp}.json
            #update taskname tag
            taskvar="rest"
            echo "$(jq --arg urlarg "${taskvar}" '. += {"TaskName": $urlarg }' ${temp}.json)" >${temp}.json
            # update skull stripped tage
            echo "$(jq '.SkullStripped = true' ${temp}.json)" >${temp}.json
            unset temp
            # clean-up
            rm -f ${headdir}/${subj}/${session}/func/denoised/func_data_smooth_usan_size.nii.gz \
                ${headdir}/${subj}/${session}/func/funcimage_BET.nii.gz \
                ${headdir}/${subj}/${session}/func/denoised/mean_func.nii.gz \
                ${headdir}/${subj}/${session}/func/denoised/report_template.html \
                ${headdir}/${subj}/${session}/func/denoised/func_data_thresh.nii.gz \
                ${headdir}/${subj}/${session}/func/denoised/mask.nii.gz
            # optional
            # rm ${headdir}/${subj}/${session}/func/denoised/${funcimagesmooth} \

            echo -e "${GREEN}--------------------------${NC}"
            echo -e "${GREEN}denoising done for ${subj}${NC}"
            echo -e "${GREEN}--------------------------${NC}"

        fi # smdesnoised exist

    fi

done
