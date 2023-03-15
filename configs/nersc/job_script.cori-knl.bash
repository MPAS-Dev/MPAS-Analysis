#!/bin/bash -l
# comment out if using debug queue
#SBATCH --partition=regular
# comment in to get premium queue
##SBATCH --qos=premium
# comment in to get the debug queue
##SBATCH --partition=debug
# comment in when run on cori haswell or knl
#SBATCH -C knl
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH -L cscratch1,SCRATCH,project

export OMP_NUM_THREADS=1

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-knl.sh
# alternatively, you can load your own development environment
# source ~/mambaforge/etc/profile.d/conda.sh
# conda activate mpas_dev
# export E3SMU_MACHINE=cori-knl

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis --polar_regions run_name_here.cfg

