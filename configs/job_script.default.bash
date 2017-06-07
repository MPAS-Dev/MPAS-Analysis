#!/bin/bash

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# no prefix is needed for jobs running on a laptop or desktop computer
command_prefix=""
# change this if not submitting this script from the directory
# containing run_analysis.py
mpas_analysis_dir="."
# the number of parallel tasks (anything between 1 and the total number
# of tasks to run)
parallel_task_count=8

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi
if [ ! -f $mpas_analysis_dir/run_analysis.py ]; then
    echo "run_analysis.py not found in $mpas_analysis_dir!"
    exit 1
fi

# This is a config file generated just for this job with the output directory,
# command prefix and parallel task count from above.
job_config_file=config.output.$RANDOM

# write out the config file specific to this job
cat <<EOF > $job_config_file
[execute]
## options related to executing parallel tasks

# the number of parallel tasks (1 means tasks run in serial, the default)
parallelTaskCount = $parallel_task_count

# Prefix on the commnd line before a parallel task (e.g. 'srun -n 1 python')
# Default is no prefix (run_analysis.py is executed directly)
commandPrefix = $command_prefix

EOF

$mpas_analysis_dir/run_analysis.py $run_config_file \
    $job_config_file

# commend this out if you want to keep the config file, e.g. for debugging
rm $job_config_file

