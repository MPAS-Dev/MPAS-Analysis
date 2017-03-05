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
  1. Create and empty config file (say `config.myrun`) or copy one of the
     example files in the `configs` directory.
  2. Copy and modify any config options you want to change from 
     `config.default` into your new config file.  Make sure they have the right
     section name (e.g. `[run]` or `[output]`).  If nothing esle, you will need
     to set `baseDirectory` under `[output]` to the folder where output should
     be stored.  Note: you should not alter `config.default` directly, since
     this is intended to hold the default configuration that is modified by your
     new config file.
  3. run: `./run_analysis.py config.myrun`.  This will read the configuraiton
     first from `config.default` and then replace that configuraiton with any
     changes from from `config.myrun`
  4. If you want to run a subset of the analysis, you can either set the 
     `generate` option under `[output]` in your config file or use the 
     `--generate` flag on the command line.  See the comments in 
     `config.default` for more details on this option.
