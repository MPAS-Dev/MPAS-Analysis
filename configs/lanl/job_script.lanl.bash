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

# change number of nodes to change the number of parallel tasks
# (anything between 1 and the total number of tasks to run)
#SBATCH --nodes=1
#SBATCH --time=1:00:00
#SBATCH --account=climateacme
#SBATCH --job-name=mpas_analysis
#SBATCH --output=mpas_analysis.o%j
#SBATCH --error=mpas_analysis.e%j
#SBATCH --qos=interactive

cd $SLURM_SUBMIT_DIR   # optional, since this is the default behavior

export OMP_NUM_THREADS=1

source /usr/projects/climate/SHARED_CLIMATE/anaconda_envs/load_latest_e3sm_unified_grizzly.sh

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
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
job_config_file=config.output.$SLURM_JOB_ID

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

mpas_analysis $run_config_file $job_config_file

