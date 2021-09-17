# MPAS-Analysis

Example config files for various HPC machines and various runs.

The intended usage is to copy or download one of these examples to the root of
MPAS-Analysis or a folder of your choice before modifying them (e.g. setting
the output `baseDirectory`) and using them to run the analysis.

To run tasks in parallel via a job queue, copy and modify the machine-specific
job script or the default job script.

## Analysis focused on polar regions

For analysis such as E3SM cryosphere simulations that is focused on polar
regions, include the `--polar_regions` flag on the command line.

Example:
```
mpas_analysis --polar_regions --generate=climatologyMapSSH --purge my_favorite.cfg 
```

