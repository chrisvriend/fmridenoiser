# fmridenoiser
 denoising pipeline for resting-state fMRI 

C. Vriend - Amsterdam UMC - Jan 9 2023

This pipeline can be used to denoise resting-state fMRI data that has been preprocessed using fmriprep. It makes use of the Denoiser tool v1.0.1 - https://github.com/arielletambini/denoiser The tool have been slightly modified to work with python 3.8 and 
requires specific versions of python packages (see requirements.txt) to run without errors, although deprecation warnings will still be produced. 

the following denoising pipelines are implemented:
24HMP8PhysSpikeReg
24HMP8PhysSpikeReg4GS
ICAAROMA8Phys
ICAAROMA8Phys4GS
for more info on these pipelines see 
https://fmridenoise.readthedocs.io/en/latest/pipelines.html
and the rsfmridenoise*.sh scripts ==> the main scripts to call.

Scripts have been optimized for use on the luna server of Amsterdam UMC (in combination with SLURM)and require several inputs (see usage info in each script).
rsfmridenoise.sh can be used to run denoising on one subject specified on the command line.
rsfmridenoise_sbatch.sh can be used to run denoising on one subject specified on the command line in combination with SLURM sbatch (sbatch ./rsfmridenoise_sbatch.sh <inputs> )
rsfmridenoise_sarray.sh can be used to run denoising on ALL subjects in the fmriprep output directory using a SLURM array.


Paths and variables will need to be changed INSIDE the script to specify the denoising pipeline, number of dummy scans to remove, smoothing kernel, etc. 




