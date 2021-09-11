#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

set -e

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate test_env
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE=compy

echo env: test_env
echo configs: ../no_polar_regions.cfg

srun -N 1 -n 1 mpas_analysis ../no_polar_regions.cfg --verbose

chmod -R ugo+rX /compyfs/www/asay932/analysis_testing/
