# MPAS-Analysis
[![Build Status](https://travis-ci.org/MPAS-Dev/MPAS-Analysis.svg?branch=develop)](https://travis-ci.org/MPAS-Dev/MPAS-Analysis)
[![Documentation Status](http://readthedocs.org/projects/mpas-analysis/badge/?version=develop)](http://mpas-analysis.readthedocs.io/en/develop/?badge=develop)


Analysis for simulations produced with Model for Prediction Across Scales
(MPAS) components and the Energy Exascale Earth System Model (E3SM), which
used those components.

![sea surface temperature](docs/_static/sst_example.png)

## Documentation

[http://mpas-analysis.readthedocs.io](http://mpas-analysis.readthedocs.io)


## Installation

MPAS-Analysis is available as an anaconda package via the `e3sm` channel:

``` bash
conda install -c conda-forge -c e3sm mpas_analysis
```

To use the latest version for developers, you will need to set up a conda
environment with the following packages:

 * numpy
 * scipy
 * matplotlib
 * netCDF4
 * xarray >= 0.10.0
 * dask
 * bottleneck
 * basemap
 * lxml
 * nco >= 4.7.0
 * pyproj
 * pillow
 * cmocean
 * progressbar2
 * requests
 * setuptools

These can be installed via the conda command:
``` bash
conda install -c conda-forge numpy scipy matplotlib netCDF4 xarray dask \
    bottleneck basemap lxml nco pyproj pillow cmocean progressbar2 requests \
    setuptools
```

## Download analysis input data

To download the data that is necessary to MPAS-Analysis, run:

``` bash
./download_analysis_data.py -o /path/to/output/directory
```

where `/path/to/output/directory` is the main folder that will contain
two subdirectories:

* `mpas_analysis`, which includes mapping and region mask files for
  standard resolution MPAS meshes
* `observations`, which includes the pre-processed observations listed in the
  [Observations table](http://mpas-analysis.readthedocs.io/en/latest/observations.html)
  and used to evaluate the model results

## List Analysis

To list the available analysis tasks, run:

``` bash
./run_mpas_analysis --list
```
This lists all tasks and their tags.  These can be used in the `generate`
command-line option or config option.  See `mpas_analysis/config.default`
for more details.

## Running the analysis

  1. Create and empty config file (say `config.myrun`), copy `config.example`,
     or copy one of the example files in the `configs` directory.
  2. Either modify config options in your new file or copy and modify config
     options from `mpas_analysis/config.default`.

  3. run: `./run_mpas_analysis config.myrun`.  This will read the configuraiton
     first from `mpas_analysis/config.default` and then replace that
     configuraiton with any changes from from `config.myrun`
  4. If you want to run a subset of the analysis, you can either set the
     `generate` option under `[output]` in your config file or use the
     `--generate` flag on the command line.  See the comments in
     `mpas_analysis/config.default` for more details on this option.

  **Requirements for custom config files:**
  * At minimum you should set `baseDirectory` under `[output]` to the folder
    where output is stored.  **NOTE** this value should be a unique
    directory for each run being analyzed.  If multiple runs are analyzed in
    the same directory, cached results from a previous analysis will not be
    updated correctly.
  * Any options you copy into the config file **must** include the
    appropriate section header (e.g. '[run]' or '[output]')
  * You do not need to copy all options from `mpas_analysis/config.default`.
    This file will automatically be used for any options you do not include
    in your custom config file.
  * You should **not** modify `mpas_analysis/config.default` directly.

## List of MPAS output files that are needed by MPAS-Analysis:

  * mpas-o files:
      * `mpaso.hist.am.timeSeriesStatsMonthly.*.nc` (Note: since OHC
        anomalies are computed wrt the first year of the simulation,
        if OHC diagnostics is activated, the analysis will need the
        first full year of `mpaso.hist.am.timeSeriesStatsMonthly.*.nc`
        files, no matter what `[timeSeries]/startYear` and
        `[timeSeries]/endYear`  are. This is especially important to know if
        short term archiving is used in the run to analyze: in that case, set
        `[input]/runSubdirectory`, `[input]/oceanHistorySubdirectory` and
        `[input]/seaIceHistorySubdirectory` to the appropriate run and archive
        directories and choose `[timeSeries]/startYear` and
        `[timeSeries]/endYear` to include only data that have been short-term
        archived).
      * `mpaso.hist.am.meridionalHeatTransport.0001-03-01.nc` (or any
        `hist.am.meridionalHeatTransport` file)
      * `mpaso.rst.0002-01-01_00000.nc` (or any other mpas-o restart file)
      * `streams.ocean`
      * `mpaso_in`
  * mpas-seaice files:
      * `mpasseaice.hist.am.timeSeriesStatsMonthly.*.nc`
      * `mpasseaice.rst.0002-01-01_00000.nc` (or any other mpas-seaice restart
        file)
      * `streams.seaice`
      * `mpassi_in`

Note: for older runs, mpas-seaice files will be named:
  * `mpascice.hist.am.timeSeriesStatsMonthly.*.nc`
  * `mpascice.rst.0002-01-01_00000.nc`
  * `streams.cice`
  * `mpas-cice_in`
  Also, for older runs mpaso-in will be named:
  * `mpas-o_in`


## Purge Old Analysis

To purge old analysis (delete the whole output directory) before running run
the analysis, add the `--purge` flag:

``` bash
./run_mpas_analysis --purge <config.file>
````
All of the subdirectories listed in `output` will be deleted along with the
climatology subdirectories in `oceanObservations` and `seaIceObservations`.

It is a good policy to use the purge flag for most changes to the config file,
for example, updating the start and/or end years of climatologies (and
sometimes time series), changing the resolution of a comparison grid, renaming
the run, changing the seasons over which climatologies are computed for a given
task, updating the code to the latest version.

Cases where it is reasonable not to purge would be, for example, changing
options that only affect plotting (color map, ticks, ranges, font sizes, etc.),
rerunning with a different set of tasks specified by the `generate` option
(though this will often cause climatologies to be re-computed with new
variables and may not save time compared with purging), generating only the
final website with `--html_only`, and re-running after the simulation has
progressed to extend time series (however, not recommended for changing the
bounds on climatologies, see above).

## Running in parallel

  1. Copy the appropriate job script file from `configs/<machine_name>` to
     the same directory as `run_mpas_analysis` (or another directory if
     preferred). The default script, `configs/job_script.default.bash`, is
     appropriate for a laptop or desktop computer with multiple cores.
  2. Modify the number of nodes (equal to the number of parallel tasks), the
     run name and optionally the output directory and the path to the config
     file for the run (default: `./configs/<machine_name>/config.<run_name>`)
     Note: in `job_script.default.bash`, the number of parallel tasks is set
     manually, since there are no nodes.
  3. Note: the number of parallel tasks can be anything between 1 and the
     number of analysis tasks to be performed.  If there are more tasks than
     parallel tasks, later tasks will simply wait until earlier tasks have
     finished.
  4. Submit the job using the modified job script

If a job script for your machine is not available, try modifying the default
job script in `configs/job_script.default.bash` or one of the job scripts for
another machine to fit your needs.


## Instructions for creating a new analysis task

Analysis tasks can be found in a directory corresponding to each component,
e.g., `mpas_analysis/ocean` for MPAS-Ocean. Shared functionality is contained
within the `mpas_analysis/shared` directory.

1. create a new task by `copying mpas_analysis/analysis_task_template.py` to
   the appropriate folder (`ocean`, `sea_ice`, etc.) and modifying it as
   described in the template.  Take a look at
   `mpas_analysis/shared/analysis_task.py` for additional guidance.
2. note, no changes need to be made to `mpas_analysis/shared/analysis_task.py`
3. modify `mpas_analysis/config.default` (and possibly any machine-specific
   config files in `configs/<machine>`)
4. import new analysis task in `mpas_analysis/<component>/__init__.py`
5. add new analysis task to `run_mpas_analysis` under `build_analysis_list`:
   ```python
      analyses.append(<component>.MyTask(config, myArg='argValue'))
   ```
   This will add a new object of the `MyTask` class to a list of analysis tasks
   created in `build_analysis_list`.  Later on in `run_analysis`, it will first
   go through the list to make sure each task needs to be generated
   (by calling `check_generate`, which is defined in `AnalysisTask`), then,
   will call `setup_and_check` on each task (to make sure the appropriate AM is
   on and files are present), and will finally call `run` on each task that is
   to be generated and is set up properly.

## Generating Documentation

To generate the `sphinx` documentation, run:
```bash
conda install sphinx sphinx_rtd_theme numpydoc m2r tabulate
cd docs
make html
```
