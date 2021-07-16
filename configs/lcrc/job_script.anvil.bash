#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH -A condo
#SBATCH -p acme-small
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh

srun -N 1 -n 1 mpas_analysis configs/polarRegions.conf 20201025.GMPAS-IAF.T62_oQU240wLI.anvil.cfg

