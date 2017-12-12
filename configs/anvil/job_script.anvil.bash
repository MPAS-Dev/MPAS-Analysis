#!/bin/bash

# comment out if using debug queue
#PBS -q acme
#PBS -A ACME
#PBS -l nodes=1
#PBS -l walltime=1:00:00
#PBS -N mpas_analysis
#PBS -o mpas_analysis.o$PBS_JOBID
#PBS -e mpas_analysis.e$PBS_JOBID

cd $PBS_O_WORKDIR

# needed to prevent interference with acme-unified
unset LD_LIBRARY_PATH
soft add +acme-unified-1.1.1-nox

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# change this if not submitting this script from the directory
# containing run_mpas_analysis
mpas_analysis_dir="."
# one parallel task per node by default
parallel_task_count=6
# ncclimo can run with 1 (serial) or 12 (bck) threads
ncclimo_mode=bck

if [ ! -f $run_config_file ]; then
    echo "File $run_config_file not found!"
    exit 1
fi
if [ ! -f $mpas_analysis_dir/run_mpas_analysis ]; then
    echo "run_mpas_analysis not found in $mpas_analysis_dir!"
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

# first, perform setup only without mpirun to create the mapping files
$mpas_analysis_dir/run_mpas_analysis --setup_only $run_config_file \
    $job_config_file
# next, do the full run now tht we have mapping files, but this time launching
# with mpirun
mpirun -n 1 $mpas_analysis_dir/run_mpas_analysis $run_config_file \
    $job_config_file

