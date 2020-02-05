#!/bin/bash
#COBALT -t 0:30:00
#COBALT -n 1
#COBALT -A OceanClimate_2
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

source /lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified.sh
export HDF5_USE_FILE_LOCKING=FALSE

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.20190301.GMPAS-DIB-IAF-ISMF.T62_oRRS30to10v3wLI.theta"

# NOTE: the following section will OVERWRITE values specified within the config file named above

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
job_config_file=config.output.$COBALT_JOBID

# write out the config file specific to this job
cat <<EOF > $job_config_file
[execute]
# options related to executing parallel tasks

# the number of parallel tasks (1 means tasks run in serial, the default)
parallelTaskCount = $parallel_task_count

# the parallelism mode in ncclimo ("serial" or "bck")
# Set this to "bck" (background parallelism) if running on a machine that can
# handle 12 simultaneous processes, one for each monthly climatology.
ncclimoParallelMode = $ncclimo_mode

EOF

# if using the mpas_analysis conda package instead of the git repo, remove
# "python -m"

python -m mpas_analysis $run_config_file $job_config_file

