#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=2:00:00
{{ sbatch }}
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

set -e

source {{ conda_base }}/etc/profile.d/conda.sh
conda activate {{ conda_env }}
export HDF5_USE_FILE_LOCKING=FALSE
export E3SMU_MACHINE={{ machine }}

echo env: {{ conda_env }}
echo configs: {{ flags }} {{ config }}

srun -N 1 -n 1 mpas_analysis --list
srun -N 1 -n 1 mpas_analysis --plot_colormaps
srun -N 1 -n 1 mpas_analysis --setup_only {{ flags }} {{ config }}
srun -N 1 -n 1 mpas_analysis --purge {{ flags }} {{ config }} --verbose
srun -N 1 -n 1 mpas_analysis --html_only {{ flags }} {{ config }}

