#!/bin/bash
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
#COBALT -t 1:00:00
#COBALT -n 1
#COBALT -A OceanClimate
#COBALT -q debug-flat-quad

export OMP_NUM_THREADS=1

module unload python
module use /projects/OceanClimate/modulefiles/all
module load e3sm-unified/1.1.2

# MPAS/ACME job to be analyzed, including paths to simulation data and
# observations. Change this name and path as needed
run_config_file="config.run_name_here"
# prefix to run a serial job on a single node on edison
command_prefix="aprun -N 1 -n 1"
# change this if not submitting this script from the directory
# containing run_mpas_analysis
mpas_analysis_dir="."
# one parallel task per node by default
parallel_task_count=12
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


# This is a config file generated just for this job with the output directory,
# command prefix and parallel task count from above.
job_config_file=config.output.$COBALT_JOBID

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

# write a batch script to prevent issues with cobalt or aprun trying to
# symlink run_mpas_analysis
cat <<EOF > run_mpas_analysis.bash
#!/bin/bash

module unload python
module use /projects/OceanClimate/modulefiles/all
module load python/anaconda-2.7-acme

$mpas_analysis_dir/run_mpas_analysis \
    $run_config_file $job_config_file
EOF

chmod u+x run_mpas_analysis.bash

$command_prefix ./run_mpas_analysis.bash

