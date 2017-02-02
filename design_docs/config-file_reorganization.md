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
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

The current example config file should be made into a general template.  Simplifications should be made to the template so that it can more easily and intuitively be modified for several analyses.  Example config files should also be added for analyzing several existing runs on several different machines.

<h2> Requirement: removal of ACME specific config options <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

To the extent possible, ACME-specific config options such as `casename` and `ref_casename_v0` should be generalized in a way that is also appropriate not just ACME runs but also any other runs involving the MPAS components we support.

<h2> Requirement: consistent section and option names <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

A consistent convention of capitalization and underscores should be used throughout the config file.


<h1> Design and Implementation </h1>

<h2> Implementation: a simple way of turning on and off individual analysis modules <br>
Date last modified: 2017/02/02 <br>
Contributors: Xylar Asay-Davis
</h2>

Implementation of the `config.template` file can be found [here](https://github.com/xylar/MPAS-Analysis/blob/5d5f64bde6ecf1d71f375a61783ff30f1654df01/config.template).


The following comment describes the planned implementation in the config file.
```
# a list of analyses to generate.  Valid names are:
#   'timeSeriesOHC', 'timeSeriesSST', 'regriddedSST', 
#   'regriddedSSS', 'regriddedMLD', 'timeSeriesSeaIceAreaVol', 
#   'regriddedSeaIceConcThick'
# the following shortcuts exist:
#   'all' -- all analyses will be run
#   'all_timeSeries' -- all time-series analyses will be run
#   'all_regriddedHorizontal' -- all analyses involving regridded horizontal
    #                                fields will be run
#   'all_ocean' -- all ocean analyses will be run
#   'all_seaIce' -- all sea-ice analyses will be run
#   'no_timeSeriesOHC' -- skip 'timeSeriesOHC' (and similarly with the
#                             other analyses).
#   'no_ocean', 'no_timeSeries', etc. -- in analogy to 'all_*', skip the 
#                                            given category of analysis
# an equivalent syntax can be used on the command line to override this
# option:
#    ./run_analysis.py config.analysis --generate \
#         all,no_ocean,all_timeSeries
generate = ['all']
```
Where there are conflicts between items in the `generate` list, successive items will override earlier items.  For example, `generate = ['all', 'no_timeSeriesOHC']` will generate all analyses except `timeSeriesOHC`.  As another example, `generate = ['all', 'no_ocean', 'all_timeSeries']` would generate all diagnostics except those comparing ocean model results with observations (and previous model results).  (Note that a more efficient and intuitive way to do the same would be `generate = ['all_seaIce', 'all_timeSeries']`.)

An analogous approach has also been added at the command line, for example:
```
./run_analysis.py config.analysis --generate all,no_ocean,all_timeSeries
```
If the `--generate` flag is used on the command line, it will replace the generate option in the config file.

As an aside, I note that it is not clear if future analysis modules will fit neatly into categories like "time series" and "regridded horizontal" fields, and these categories are not meant to be all-encompassing.

<h2> Implementation: there should be a simplified template for config files <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

The required `config.template` has been implemented in #86, specifically [here](https://github.com/xylar/MPAS-Analysis/blob/5d5f64bde6ecf1d71f375a61783ff30f1654df01/config.template).  A subdirectory `configs` will be added with several examples from runs on LANL IC and on Edison at NERSC.  Other examples can be added as appropriate and useful.

<h2> Implementation: removal of ACME specific config options <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

`casename` has been renamed `mainRunName`, `referenceRunName` has been added for comparison with reference runs that have not been preprocessed (not yet supported), and `ref_casename_v0` has been renamed `preprocessedReferenceRunName`.

See #86, specifically [config.template](https://github.com/xylar/MPAS-Analysis/blob/5d5f64bde6ecf1d71f375a61783ff30f1654df01/config.template).

<h2> Implementation: consistent section and option names <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

In [config.template](https://github.com/xylar/MPAS-Analysis/blob/5d5f64bde6ecf1d71f375a61783ff30f1654df01/config.template) in #86, "[CamelCase](https://en.wikipedia.org/wiki/Camel_case)" has been used for all sections and options.  The first word is lowercase and subsequent words begin with an uppercase latter.  Underscores have been removed (except in the syntax used to turn on and off options, where underscores in prefixes `all_` and `no_` make splitting and comparison simpler in the implementation.


<h1> Testing </h1>

<h2> Testing and Validation: a simple way of turning on and off individual analysis modules <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

CI will be added to make sure that the function to parse the generate list (`run_analysis.check_generate`) behaves as expected.

<h2> Testing and Validation: there should be a simplified template for config files <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

There is not a way to test the template in the usual sense.  Instead, the test will be asking other developers and users to adapt the template for new runs to make sure it is intuitive.

<h2> Testing and Validation: removal of ACME specific config options <br>
Date last modified: 2017/01/29 <br>
Contributors: Xylar Asay-Davis
</h2>

For now, the plan is just to rename the appropriate config options, so the test is simply to ensure that analysis runs correctly and produces bit-for-bit identical images to those produced by the current `MPAS-Analysis/develop`.


<h2> Testing and Validation: consistent section and option names <br>
Date last modified: 2017/02/01 <br>
Contributors: Xylar Asay-Davis
</h2>

As above, the test is simply to ensure that analysis runs correctly and produces bit-for-bit identical images to those produced by the current `MPAS-Analysis/develop`.


