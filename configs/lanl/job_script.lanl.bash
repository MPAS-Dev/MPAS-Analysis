#!/bin/bash
#SBATCH  --job-name=mpas_analysis
#SBATCH  --nodes=1
#SBATCH  --output=mpas_analysis.o%j
#SBATCH  --time=1:00:00
#SBATCH  --qos=standard
#SBATCH  --partition=standard

source /users/xylar/climate/mambaforge/etc/profile.d/conda.sh
source /users/xylar/climate/mambaforge/etc/profile.d/mamba.sh
mamba activate mpas_dev

export HDF5_USE_FILE_LOCKING=FALSE

mpas_analysis --verbose -m chicoma-cpu --polar_regions 20230109.GMPAS-JRA1p5-DIB-ISMF.TL319_SOwISC12to60E2r4.chicoma-cpu.cfg --purge
