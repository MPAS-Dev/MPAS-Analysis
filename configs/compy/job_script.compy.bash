#!/bin/bash -l
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=e3sm
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1

source /compyfs/software/e3sm-unified/load_latest_e3sm_unified.sh
export HDF5_USE_FILE_LOCKING=FALSE

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi

# if using the mpas_analysis conda package instead of the git repo, remove
# "python -m"

srun -N 1 -n 1 python -m mpas_analysis $run_config_file

