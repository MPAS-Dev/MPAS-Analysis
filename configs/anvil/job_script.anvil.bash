#!/bin/bash
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

# comment out if using debug queue
#PBS -q acme
#PBS -A ACME
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -N mpas_analysis
#PBS -o mpas_analysis.o$PBS_JOBID
#PBS -e mpas_analysis.e$PBS_JOBID

cd $PBS_O_WORKDIR

source /lcrc/soft/climate/e3sm-unified/base/etc/profile.d/conda.sh
conda activate e3sm_unified_1.2.0_py2.7_nox
export HDF5_USE_FILE_LOCKING=FALSE
# needed to prevent interference with acme-unified
unset LD_LIBRARY_PATH


# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# number of parallel tasks to run
parallel_task_count=6
# ncclimo can run with 1 (serial) or 12 (bck) threads
ncclimo_mode=bck

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi

job_config_file=config.output.$PBS_JOBID

# write out the config file specific to this job
cat <<EOF > $job_config_file
[execute]
## options related to executing parallel tasks

# the number of parallel tasks (1 means tasks run in serial, the default)
parallelTaskCount = $parallel_task_count

# the parallelism mode in ncclimo ("serial" or "bck")
# Set this to "bck" (background parallelism) if running on a machine that can
# handle 12 simultaneous processes, one for each monthly climatology.
ncclimoParallelMode = $ncclimo_mode

EOF

# Note: to run from the conda package rather than a local git repo, remove
# "python -m" in each command below

# first, perform setup only without mpirun to create the mapping files
python -m mpas_analysis --setup_only $run_config_file $job_config_file
# next, do the full run now that we have mapping files, but this time launching
# with mpirun
mpirun -n 1 python -m mpas_analysis $run_config_file $job_config_file
