#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH -A cli115
#SBATCH -p batch
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

source /gpfs/alpine/proj-shared/cli115/e3sm-unified/load_latest_e3sm_unified_andes.csh
# alternatively, you can load your own development environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate mpas_dev
# export E3SMU_MACHINE=anvil

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis run_name_here.cfg

