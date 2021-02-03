#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /home/ac.xylar/miniconda3/etc/profile.d/conda.sh
conda activate test_env
export HDF5_USE_FILE_LOCKING=FALSE

echo env: test_env
echo configs: ../../configs/polarRegions.conf ../main.cfg

mpas_analysis --list
mpas_analysis --plot_colormaps
mpas_analysis --setup_only ../../configs/polarRegions.conf ../main.cfg
mpas_analysis --purge ../../configs/polarRegions.conf ../main.cfg --verbose
mpas_analysis --html_only ../../configs/polarRegions.conf ../main.cfg

