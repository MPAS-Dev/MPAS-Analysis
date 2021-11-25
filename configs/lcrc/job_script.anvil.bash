#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH -A condo
#SBATCH -p acme-small
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

export OMP_NUM_THREADS=1

source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh
# alternatively, you can load your own development environment
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate mpas_dev
# export E3SMU_MACHINE=anvil

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis --polar_regions 20201025.GMPAS-IAF.T62_oQU240wLI.anvil.cfg

