#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

export OMP_NUM_THREADS=1

source ~/mambaforge/etc/profile.d/conda.sh
conda activate mpas_analysis_dev
# if you are on an E3SM supported machine, you can specify it:
# export E3SMU_MACHINE=chrysalis

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis --polar_regions run_name_here.cfg
