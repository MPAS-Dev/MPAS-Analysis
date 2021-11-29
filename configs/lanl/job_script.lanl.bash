#!/bin/bash
# change number of nodes to change the number of parallel tasks
# (anything between 1 and the total number of tasks to run)
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH --qos=interactive

export OMP_NUM_THREADS=1

source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/load_latest_e3sm_unified_grizzly.sh
# alternatively, you can load your own development environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate mpas_dev
# export E3SMU_MACHINE=grizzly

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis run_name_here.cfg

