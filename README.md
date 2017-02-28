# MPAS-Analysis
[![Build Status](https://travis-ci.org/MPAS-Dev/MPAS-Analysis.svg?branch=develop)](https://travis-ci.org/MPAS-Dev/MPAS-Analysis)

A repository for the development and maintenance of MPAS analysis tools.

Analysis is stored in a directory corresponding to each core component, e.g., `ocean` for
MPAS-Ocean. Shared functionality is contained within the `shared` directory.

## Installation
This analysis repository presumes that the following python packages are available:

 * numpy
 * scipy
 * matplotlib
 * netCDF4
 * xarray
 * dask
 * bottleneck
 * basemap
 * nco with `conda install -c conda-forge nco`

You can easily install them via the conda command:

```
conda install -c conda-forge numpy scipy matplotlib netCDF4 xarray dask bottleneck basemap nco
```

## Running the analysis
  1. create a configuration file by copying and modifying `config.template` or one of the example files in `configs`
  2. run: `./run_analysis.py config.myrun`, where `config.myrun` is the config file you created
  3. If you want to run a subset of the analysis, you can either modify the `generate` option under `[output]` in the config file or use the `--generate` flag on the command line.  See `config.template` for more details.
