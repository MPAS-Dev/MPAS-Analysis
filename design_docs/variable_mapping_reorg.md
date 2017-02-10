<h1> Title: Moving variable mapping outside of mpas_xarray <br>
Xylar Asay-Davis <br>
date: 2017/02/10 <br>
</h1>
<h2> Summary </h2>
In discussions with @pwolfram, it became clear that we would like to keep
mpas_xarray as general as possible, rather than adding code specific to
MPAS-Analysis.  In particular, the capability for mapping variable names
that is currently part of mpas_xarray is likely a capability that only
MPAS-Analysis will need when opening xarray data sets.  Likewise, there is
a desire for mpax_xarray not to use any of the functionality outside of its
own module so that it remains autonomous from MPAS-Analysis.

At the same time, it is desirable for efficiency and parallelism to perform
certain operations during the preprocessing step within xarray, rather than
constructing a data set first and then (in serial) performing manipulations
(e.g. creating a time coordinate and slicing variables).

The solution will be tested by making sure it produces bit-for-bit identical
results to those from the develop branch for typical test cases on LANL IC
and Edison.

<h1> Requirements </h1>

<h2> Requirement: mpas_xarray does not include MPAS-Analysis specific
functionality <br>
Date last modified: 2017/02/10 <br>
Contributors: Xylar Asay-Davis
</h2>

MPAS-Analysis specific functionality such as variable mapping should be 
removed from mpas_xarray so it can remain an independent module, requiring
minimal modification to accommodate MPAS-Analysis' needs.

<h2> Requirement: MPAS-Analysis specific functionality should be supported in
xarray preprossing <br>
Date last modified: 2017/02/10 <br>
Contributors: Xylar Asay-Davis
</h2>

There should be a way to perform MPAS-Analysis specific functionality such as
mapping variables during preprocessing.  This functionality should be
relatively easy to add to as new preprocessing needs arise.


<h1> Algorithmic Formulations (optional) </h1>

<h2> Algorithm: mpas_xarray does not include MPAS-Analysis specific
functionality <br>
Date last modified: 2017/02/10 <br>
Contributors: Xylar Asay-Davis
</h2>

All functions and function arguments related to variable mapping will
be removed from mpas_xarray and moved elsewhere.

<h2> Algorithm: MPAS-Analysis specific functionality should be supported in
xarray preprossing <br>
Date last modified: 2017/02/15 <br>
Contributors: Xylar Asay-Davis
</h2>

A new utility function, `open_multifile_dataset` will added to `mpas_xarray` 
that simplifies current calls to `xarray.open_mfdataset` to hide the 
preprocessor and take care of removing redundant time indices once the dataset
has been built. (This function doesn't directly address the requirement but
is meant to make `mpas_xarray` easier to use and made sense because it
has a one-to-one correspondence with other functionality, described below,
that does address the requirement.)

A new module, `generalized_reader` will also be added with its own
`open_multifile_dataset` function.  This version takes additional arguments
including a variable map and start and end dates for the dataset.
`generalized_reader.open_multifile_dataset` will create a data set
by calling `xarray.open_mfdataset` with its own preprocessing function,
`generalized_reader._preprocess` that first maps variable names, then
calls `mpas_xarray.preprocess` to finish the job.  Once the dataset has 
been constructed, redundant time indices are removed and the 'Time'
coordinate is sliced to be between the supplied start and end dates.

This solution may add some confusion in terms of which reader should
be used to open xarray datasets.  It is my sense that most developers
adding new functionality will do so by modifying existing scripts, and
these examples should make it clear which version of 
`open_multifile_dataset` is most appropriate.  Nevertheless, clear
documentation of `generalized_reader` and `mpas_xarray`, and their
differences are needed.

Here is a typical usage of `generalized_reader.open_multifile_dataset`:
```python
from mpas_analysis.shared.generalized_reader.generalized_reader \
    import open_multifile_dataset
    
file_name = 'example_jan_feb.nc'
timestr = ['xtime_start', 'xtime_end']
var_list = ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']
variable_map = {
   'avgSurfaceTemperature':
       ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature',
        'other_string',
        'yet_another_string'],
   'daysSinceStartOfSim':
       ['time_avg_daysSinceStartOfSim',
        'xtime',
        'something_else']}
ds = open_multifile_dataset(file_names=file_name,
                            calendar=calendar,
                            time_variable_name=timestr,
                            variable_list=var_list,
                            start_date='0001-01-01',
                            end_date='9999-12-31',
                            variable_map=variable_map,
                            year_offset=1850)
```

Here is the same for `mpas_xarray.open_multifile_dataset` without the 
variable map, start and end dates:
```python
from mpas_analysis.shared.mpas_xarray.mpas_xarray \
    import open_multifile_dataset
    
file_name = 'example_jan_feb.nc'
timestr = ['xtime_start', 'xtime_end']
var_list = ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

ds = open_multifile_dataset(file_names=file_name,
                            calendar=calendar,
                            time_variable_name=timestr,
                            variable_list=var_list,
                            year_offset=1850)
```


<h1> Design and Implementation </h1>

<h2> Implementation: mpas_xarray does not include MPAS-Analysis specific
functionality <br>
Date last modified: 2017/02/15 <br>
Contributors: Xylar Asay-Davis
</h2>

A test branch can be found here 
[xylar/MPAS-Analysis/variable_mapping_reorg](https://github.com/xylar/MPAS-Analysis/tree/variable_mapping_reorg)

I have removed `map_variable` and `rename_variables` from `mpas_xarray`.
I also removed any mention of the variable map from the rest of `mpas_xarray`.

This branch also includes several other cleanup operations that are not
addressing any requirements.  These include:
   - I added a new helper function, `open_multifile_dataset`, for opening an
     xarray data set in a single, simple command without reference to the
     preprocessor.  This function should make opening new data sets more
     intuitive for mpas_xarray users.
   - making several utility functions non-public (it is unclear to me why anyone
     want to call these directly):
       - `_assert_valid_datetimes`
       - `_assert_valid_selections`
       - `_ensure_list`
       - `_get_datetimes`
   - I have removed the ability to run `mpas_xarray.py` as a script and the associated
     tests.  This is on the premise that 1) the test were outdated and would have
     needed to be updated to work with the current code and 2) unit testing in
     `test/test_mpas_xarray.py` takes care of this capability in a better way.
   - I have tried to make variable names a bit more verbose in various places.
     However, at @pwolfram'2 request, I have left ds for datasets, following the
     `xarray` convention.
   - I have tried to improve the docstrings using a syntax that should be useful
     for generating documentation later on.
   - I have update unit testing to work with the new inerface, notably the
     `open_multifile_dataset` function.

<h2> Implementation: MPAS-Analysis specific functionality should be supported in
xarray preprossing <br>
Date last modified: 2017/02/15 <br>
Contributors: Xylar Asay-Davis
</h2>

In the same branch as above, I have added a `generalized_reader` module that
extends the capabilities of `mpas_xarray` to include mapping of variable names.
The file structure is as follows:
```
mpas_analysis/shared/
             -  generalized_reader/
                     __init__.py
                    generalized_reader.py
```
`generalized_reader.py` contains a function `open_multifile_dataset` that is similar to
the one in `mpas_xarray` but with additional arguments needed by analysis:
  - `variable_map`, a map between MPAS and MPAS-Analysis variable names
  - `start_date`, the start date of the analysis
  - `end_date`, the end date of the analysis
This function performs the same steps as `mpas_xarray.open_multifile_dataset`
but uses the local preprocessing function, `_preprocess`, and also slices
the 'Time' coordinate using the given start and end dates as a final step.

The `generalized_reader._preprocess` funciton first maps variable names, then calls 
`mpas_xarray.preprocess` to do the rest of the preprocessing as normal.

Two private functions, `_map_variable_name` and `_rename_variables` (take out of
`mpas_xarray`) are used to perform variable-name mapping.

<h1> Testing </h1>

<h2> Testing and Validation: MPAS-Analysis specific functionality should be supported in
xarray preprossing <br>
Date last modified: 2017/02/15 <br>
Contributors: Xylar Asay-Davis
</h2>

In [xylar/MPAS-Analysis/variable_mapping_reorg](https://github.com/xylar/MPAS-Analysis/tree/variable_mapping_reorg),
the unit testing for mpas_xarray has been updated.  This includes moving unit testing for
variable mapping elsewhere.

I will make sure all tests with config files in the `configs/lanl` and `configs/edison` 
directories produce bit-for-bit results with the current `develop`.

<h2> Testing and Validation: MPAS-Analysis specific functionality should be supported in
xarray preprossing <br>
Date last modified: 2017/02/10 <br>
Contributors: Xylar Asay-Davis
</h2>

Largely, the same as above.

I have added unit testing for `generalized_reader` (via the standalone
`generalized_reader.open_multifile_dataset` function).  These tests ensure that:
  - variable mapping works as expected
  - start and end dates work as expected
