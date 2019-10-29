#!/bin/bash -l
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

# comment out if using debug queue
#SBATCH --partition=regular
# comment in to get premium queue
##SBATCH --qos=premium
# comment in to get the debug queue
##SBATCH --partition=debug
# comment in when run on cori haswell or knl
#SBATCH -C knl
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=acme
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH -L cscratch1,SCRATCH,project

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

export OMP_NUM_THREADS=1

source /global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified.sh
export HDF5_USE_FILE_LOCKING=FALSE

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi

# For an E3SM cryosphere run, include configs/polarRegions.conf, or exclude
# this extra config file for defalut parameters
srun -N 1 -n 1 python -m mpas_analysis configs/polarRegions.conf $run_config_file

# if using the mpas_analysis conda package instead of the git repo
# srun -N 1 -n 1 mpas_analysis $run_config_file

