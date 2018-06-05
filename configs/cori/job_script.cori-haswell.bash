#!/bin/bash -l
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

# comment out if using debug queue
#SBATCH --partition=regular
# comment in to get premium queue
##SBATCH --qos=premium
# comment in to get the debug queue
##SBATCH --partition=debug
# comment in when run on cori haswell or knl
#SBATCH -C haswell
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=acme
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH -L cscratch1,SCRATCH,project

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

export OMP_NUM_THREADS=1

module unload python python/base
module use /global/project/projectdirs/acme/software/modulefiles/all
module load e3sm-unified/1.1.2

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi
if [ ! -f ./run_mpas_analysis ]; then
    echo "run_mpas_analysis not found in current directory!"
    exit 1
fi

srun -N 1 -n 1 ./run_mpas_analysis $run_config_file

