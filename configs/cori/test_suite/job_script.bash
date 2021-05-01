#!/bin/bash -l
#SBATCH --partition=regular
#SBATCH -C haswell
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH -L cscratch1,SCRATCH,project

set -e

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate test_env
export HDF5_USE_FILE_LOCKING=FALSE

echo env: test_env
echo configs: ../../configs/polarRegions.conf ../main.cfg

srun -N 1 -n 1 mpas_analysis --list
srun -N 1 -n 1 mpas_analysis --plot_colormaps
srun -N 1 -n 1 mpas_analysis --setup_only ../../configs/polarRegions.conf ../main.cfg
srun -N 1 -n 1 mpas_analysis --purge ../../configs/polarRegions.conf ../main.cfg --verbose
srun -N 1 -n 1 mpas_analysis --html_only ../../configs/polarRegions.conf ../main.cfg

chmod -R ugo+rX /global/cfs/cdirs/e3sm/www/xylar/analysis_testing/
