#!/bin/bash
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

# comment out if using debug queue
#PBS -q batch
# comment in to get the debug queue (only available on Titan)
##PBS -q debug
# change number of nodes to change the number of parallel tasks
# (anything between 1 and the total number of tasks to run)
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -A cli115
#PBS -N mpas_analysis
#PBS -o mpas_analysis.o$PBS_JOBID
#PBS -e mpas_analysis.e$PBS_JOBID

cd $PBS_O_WORKDIR

source /ccs/proj/cli900/sw/rhea/e3sm-unified/load_latest_e3sm_unified.csh
export HDF5_USE_FILE_LOCKING=FALSE

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# command to run a serial job on a single node on olcf machines.  Remove
# "python -m" is using the conda package instead of a git repo
command="aprun -b -N 1 -n 1 python -m mpas_analysis"
# one parallel task per node by default
parallel_task_count=12
# ncclimo can run with 1 (serial) or 12 (bck) threads
ncclimo_mode=bck

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi

# This is a config file generated just for this job with the output directory,
# command prefix and parallel task count from above.
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

$command $run_config_file $job_config_file

