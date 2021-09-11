#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

set -e

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /home/ac.xylar/chrysalis/miniconda3/etc/profile.d/conda.sh
conda activate test_env
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE=chrysalis

echo env: test_env
echo configs: ../../configs/polar_regions.cfg ../main.cfg

srun -N 1 -n 1 mpas_analysis --list
srun -N 1 -n 1 mpas_analysis --plot_colormaps
srun -N 1 -n 1 mpas_analysis --setup_only ../../configs/polar_regions.cfg ../main.cfg
srun -N 1 -n 1 mpas_analysis --purge ../../configs/polar_regions.cfg ../main.cfg --verbose
srun -N 1 -n 1 mpas_analysis --html_only ../../configs/polar_regions.cfg ../main.cfg

