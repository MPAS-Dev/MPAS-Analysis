#!/bin/bash

# comment out if using debug queue
#PBS -q batch
# comment in to get the debug queue (only available on Titan)
##PBS -q debug
# change number of nodes to change the number of parallel tasks
# (anything between 1 and the total number of tasks to run)
#PBS -l nodes=10
#PBS -l walltime=1:00:00
#PBS -A cli115
#PBS -N mpas_analysis
#PBS -o mpas_analysis.o$PBS_JOBID
#PBS -e mpas_analysis.e$PBS_JOBID

cd $PBS_O_WORKDIR

module unload python
module use /ccs/proj/cli115/pwolfram/modulefiles/all
module load python/anaconda-2.7-climate

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# prefix to run a serial job on a single node on edison
command_prefix="aprun -b -N 1 -n 1"
# change this if not submitting this script from the directory
# containing run_mpas_analysis
mpas_analysis_dir="."
# one parallel task per node by default
parallel_task_count=$PBS_NUM_NODES

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi
if [ ! -f $mpas_analysis_dir/run_mpas_analysis ]; then
    echo "run_mpas_analysis not found in $mpas_analysis_dir!"
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

# Prefix on the commnd line before a parallel task (e.g. 'srun -n 1 python')
# Default is no prefix (run_mpas_analysis is executed directly)
commandPrefix = $command_prefix

EOF

$mpas_analysis_dir/run_mpas_analysis $run_config_file \
    $job_config_file

