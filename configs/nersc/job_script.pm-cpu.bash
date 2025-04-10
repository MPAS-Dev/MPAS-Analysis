#!/bin/bash
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH --account=e3sm
#SBATCH --nodes=2
#SBATCH --exclusive
#SBATCH --time=1:00:00
#SBATCH --qos=regular
#SBATCH --constraint=cpu

export OMP_NUM_THREADS=1

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
# alternatively, you can load your own development environment
# source ~/mambaforge/etc/profile.d/conda.sh
# conda activate mpas_analysis_dev
# export E3SMU_MACHINE=pm-cpu

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis --polar_regions run_name_here.cfg
