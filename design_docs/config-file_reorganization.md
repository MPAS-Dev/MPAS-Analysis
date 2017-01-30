<h1> Config File Reorganization <br>
Xylar Asay-Davis<br>
date: 01-29-2017<br>
</h1>
<h2> Summary </h2>
This document describes various efforts to clean up the structure of the MPAS-Analysis config file.  The idea is to create a template config file that will replace `config.analysis` as well as a number of example config files designed to make use of various MPAS and ACME runs on various machines.  The reorganization should make the analysis easier for users to modify and run.

<h1> Requirements </h1>

<h2> Requirement: a simple way of turning on and off individual analysis modules <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

There should be a simple, intuitive method for turning on and off individual analysis modules (e.g. `ocean/ohc_timeseries`).  This should replace the current approach of having a boolean `generate` flag for each analysis module in a separate config section.  Preferably, there should be an equivalent method for turning on and off analysis modules from the command line that overrides that in the config file.

<h2> Requirement: there should be a simplified template for config files <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

The current example config file is specific to a run on Edison and should be made into a general template.  Simplifications should be made to the template so that it can more easily and intuitively be modified for several analyses.  Example config files should also be added for analyzing several existing runs on several different machines.

<h2> Requirement: removal of ACME specific config options <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

To the extent possible, ACME-specific config options such as `casename` and `ref_casename_v0` should be generalized in a way that is also appropriate for MPAS-only runs.

<h2> Requirement: optional reference dates for each model run or set of observations<br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

Currently, there is an option `yr_offset` that is intended to apply to all dates and observations.  This option should be removed and possibly be replaced by an input and output reference date for each model run and observational data set.  Discussion is needed if these reference dates are actually needed once we find a more generalized calendar (e.g. `netcdftime`).

<h1> Design and Implementation </h1>

<h2> Implementation: a simple way of turning on and off individual analysis modules <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

The following comment describes the planned implementation in the config file.
```
# a list of which analyses to generate.  Valid names are:
#   'ohc_timeseries', 'sst_timeseries', 'sst_modelvsobs', 'sss_modelvsobs',
#   'mld_modelvsobs', 'seaice_timeseries', 'seaice_modelvsobs'
# the following shortcuts exist:
#   'all' -- all analyses will be run
#   'all_timeseries' -- all time-series analyses will be run
#   'all_modelvsobs' -- all analyses comparing model climatologies with
#                       observations will be run
#   'all_ocean' -- all ocean analyses will be run
#   'all_seaice' -- all sea-ice analyses will be run
#   'no_ohc_timeseries' -- skip the 'ohc_timeseries' (and similarly with the other analyses).
#   'no_ocean', 'no_timeseries', etc. -- in analogy to 'all_*', skip the given category of analysis
generate = ['all']
```
Where there are conflicts between items in the `generate` list, successive items will override earlier items.  For example, `generate = ['all', 'no_ohc_timeseries']` will generate all analyses except `ohc_timeseries`.  As another example, `generate = ['all', 'no_ocean', 'all_timeseries']` would generate all diagnostics except those comparing ocean model results with observations (and previous model results).  (Note that a more efficient and intuitive way to do the same would be `generate = ['all_seaice', 'all_timeseries']`.)

An analogous approach will also be added at the command line, for example:
```
./run_analysis.py config.analysis --generate all,no_ocean,all_timeseries
```
(I am open to other syntax.)  If the --generate flag is used on the command line, it will replace the generate option in the config file.

As an aside, I note that it is not clear if future analysis modules will fit neatly into categories like "time series" and "model vs. observations", and these categories are not meant to be all-encompassing.

<h2> Implementation: there should be a simplified template for config files <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

Such a template has been implemented in #86.  A subdirectory `configs` will be added with several examples from runs on LANL IC and on Edison at NERSC.  Other examples can be added as appropriate and useful.

<h2> Implementation: removal of ACME specific config options <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

This item needs some discussion.  in #86, I would like to renamed `ref_casename_v0` to simply `ref_casename`.  `casename` and `ref_casename_v0` are used for file names and figure titles (essentially legends).  It would be useful to discuss what the relevant equivalents would be fore standalone MPAS runs.

<h1> Testing </h1>

<h2> Testing and Validation: a simple way of turning on and off individual analysis modules <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

This will be difficult to test with CI since we currently don't have CI that can perform tests with `run_analysis.py`.

Instead, I will manually test a variety of combinations of `generate` lists (both in the config file and on the command line).  I will list the tests I perform in #86.

<h2> Testing and Validation: there should be a simplified template for config files <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

There is not a way to test the template in the usual sense.  Instead, the test will be asking other developers and users to adapt the template for new runs to make sure it is intuitive.

<h2> Testing and Validation: removal of ACME specific config options <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

The analysis will be tested on both MPAS and ACME runs to make sure the resulting image files have useful file names and titles (legends).  
