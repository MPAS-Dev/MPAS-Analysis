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

source /home/xylar/miniconda3/etc/profile.d/conda.sh
conda activate mpas-analysis
export HDF5_USE_FILE_LOCKING=FALSE

srun -N 1 -n 1 python -m mpas_analysis configs/polarRegions.conf config.SO60to10ISC.20200108

