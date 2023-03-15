#!/bin/bash
#COBALT -t 3:00:00
#COBALT -n 1
#COBALT -A ClimateEnergy_4

source /lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified_cooley.sh
# alternatively, you can load your own development environment
# source ~/mambaforge/etc/profile.d/conda.sh
# conda activate mpas_dev
# export E3SMU_MACHINE=cooley

export HDF5_USE_FILE_LOCKING=FALSE

# For an E3SM cryosphere run, include --polar_regions, or exclude
# this extra flag for default parameters
mpas_analysis 20180410.A_WCYCL1950_HR.ne120_oRRS18v3_ICG.theta.cfg

