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
 * xarray ≥ 0.9.1
 * dask
 * bottleneck
 * basemap
 * lxml
 * nco

You can easily install them via the conda command:

```
conda config --add channels conda-forge
conda install numpy scipy matplotlib netCDF4 xarray dask bottleneck basemap lxml nco
```

## Running the analysis
  1. Create and empty config file (say `config.myrun`) or copy one of the
     example files in the `configs` directory.
  2. Copy and modify any config options you want to change from
     `config.default` into your new config file.

     **Requirements for custom config files:**
     * At minimum you should set `baseDirectory` under `[output]` to the folder
       where output is stored.  **NOTE** this value should be a unique
       directory for each run being analyzed.  If multiple runs are analyzed in
       the same directory, cached results from a previous analysis will not be
       updated correctly.
     * Any options you copy into the config file **must** include the
       appropriate section header (e.g. '[run]' or '[output]')
     * The entire `config.default` does not need to be used.  This fill will
       automatically be used for any options you do not include in your custom
       config file.
     * Given the automatic sourcing of `config.default` you should **not**
       alter `config.default` directly.
  3. run: `./run_analysis.py config.myrun`.  This will read the configuraiton
     first from `config.default` and then replace that configuraiton with any
     changes from from `config.myrun`
  4. If you want to run a subset of the analysis, you can either set the
     `generate` option under `[output]` in your config file or use the
     `--generate` flag on the command line.  See the comments in
     `config.default` for more details on this option.
