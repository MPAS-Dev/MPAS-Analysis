#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

srun --mpi=pmi2 -N 1 -n 1 mpas_analysis config.run_name_here

