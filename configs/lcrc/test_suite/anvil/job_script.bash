#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH -A condo
#SBATCH -p acme-small
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

set -e

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /home/ac.xylar/anvil/miniconda3/etc/profile.d/conda.sh
conda activate test_env
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE=anvil

echo env: test_env
echo configs: --polar_regions ../main.cfg

srun -N 1 -n 1 mpas_analysis --list
srun -N 1 -n 1 mpas_analysis --plot_colormaps
srun -N 1 -n 1 mpas_analysis --setup_only --polar_regions ../main.cfg
srun -N 1 -n 1 mpas_analysis --purge --polar_regions ../main.cfg --verbose
srun -N 1 -n 1 mpas_analysis --html_only --polar_regions ../main.cfg

