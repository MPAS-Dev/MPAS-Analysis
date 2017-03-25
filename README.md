# MPAS-Analysis
[![Build Status](https://travis-ci.org/MPAS-Dev/MPAS-Analysis.svg?branch=develop)](https://travis-ci.org/MPAS-Dev/MPAS-Analysis)

A repository for the development and maintenance of MPAS analysis tools.

![sea surface temperature](docs/_static/sst_example.png)

Analysis is stored in a directory corresponding to each core component, e.g., `ocean` for
MPAS-Ocean. Shared functionality is contained within the `shared` directory.

## Installation
This analysis repository presumes that the following python packages are available:

 * numpy
 * scipy
 * matplotlib
 * netCDF4
 * xarray >= 0.9.1
 * dask
 * bottleneck
 * basemap
 * lxml
 * nco >= 4.6.8
 * pyproj
 * pillow

You can easily install them via the conda command:

```
conda config --add channels conda-forge
conda install numpy scipy matplotlib netCDF4 xarray dask bottleneck basemap \
    lxml nco pyproj
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


## Running in parallel
  1. Copy the appropriate job script file from `configs/<machine_name>` to
     the same directory as `run_analysis.py` (or another directory if preferred).
     The default script, `configs/job_script.default.bash`, is appropriate for
     a laptop or desktop computer with multiple cores.
  2. Modify the number of nodes (equal to the number of parallel tasks), the
     run name and optionally the output directory and the path to the config
     file for the run (default: `./configs/<machine_name>/config.<run_name>`)
     Note: in `job_script.default.bash`, the number of parallel tasks is set
     manually, since there are no nodes.
  3. Note: the number of parallel tasks can be anything between 1 and the number
     of analysis tasks to be performed.  If there are more tasks than parallel
     tasks, later tasks will simply wait until earlier tasks have finished.
  4. Submit the job using the modified job script

If a job script for your machine is not available, try modifying the default
job script in `configs/job_script.default.bash` or one of the job scripts for
another machine to fit your needs.


## Instructions for creating a new analysis task

1. create a new task by `copying mpas_analysis/analysis_task_template.py` to
   the appropriate folder (`ocean`, `sea_ice`, etc.) and modifying it as
   described in the template.  Take a look at
   `mpas_analysis/shared/analysis_task.py` for additional guidance.
2. note, no changes need to be made to `mpas_analysis/shared/analysis_task.py`
3. modify `config.default` (and possibly any machine-specific config files in
   `configs/<machine>`)
4. import new analysis task in `mpas_analysis/<component>/__init__.py`
5. add new analysis task to `run_analysis.py` under `build_analysis_list`:
   ```python
      analyses.append(<component>.MyTask(config, myArg='argValue'))
   ```
   This will add a new object of the `MyTask` class to a list of analysis tasks
   created in `build_analysis_list`.  Later on in `run_analysis`, it will first
   go through the list to make sure each task needs to be generated
   (by calling `check_generate`, which is defined in `AnalysisTask`), then, will
   call `setup_and_check` on each task (to make sure the appropriate AM is on
   and files are present), and will finally call `run` on each task that is
   to be generated and is set up properly.

## Generating Documentation

To generate the `sphinx` documentation, run:
```bash
conda install sphinx sphinx_rtd_theme numpydoc
cd docs
make html
```
