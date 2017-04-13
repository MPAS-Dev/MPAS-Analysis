<h1> Title: Generalize Calendar supported by Analysis <br>
Xylar Asay-Davis <br>
date: 2017/02/09 <br>
</h1>
<h2> Summary </h2>
Currently, the time variable in `xarray` data sets within MPAS-Analysis has two
major shortcomings, inherited from `xarray` (through `pandas` and `numpy.datetime64`).
First, only the Gregorian calendar is supported.  Second, there is not support
for dates outside the years 1678 to 2262.  The analysis needs to support both
the Gregorian  ('gregorian') and the 365-day ('gregorian_noleap') calendars.  It also needs to
support, at a minimum, years between 0001 and 9999, and preferably arbitrary
years both positive and negative.

A major challenge is that it seems that xarray cannot easily be forced to
use an alternative representation of dates to the troublesome
`numpy.datetime64` type (see, for example, 
[pydata/xarray#1084](https://github.com/pydata/xarray/issues/1084)).  
The most obvious alternative, `datetime.datetime`,
seemingly cannot be used directly in `xarray` because objects of this type
are converted to `numpy.datetime64` objects at various stages when using 
features from pandas, raising errors when dates are out of range.  While an
alternative date class (e.g. `netcdftime.DatetimNoLeap`) might be used to
represent dates on the 'gregorian_noleap' calendar, there is no such
preexisting alternative for the 'gregorian' calendar.

The solution proposed herein is to store time as floating-point days since the
reference date 0001-01-01 and to convert dates in this format to 
`datetime.datetime` and `MpasRelativeDelta` objects whenever mathematical
manipulation of dates is required.

A successful implementation would produce essentially identical analysis to
what is currently produced, but making use of the dates from the MPAS calendar
(whether Gregorian or 365-day) without the need for artifical offsets (e.g.
`yearOffset` used in the current code.  Plots of horizontal fields would remain
unchanged while plots of time series would have a time axis with the simulation
date instead of the offset date.


<h1> Requirements </h1>

<h2> Requirement: The 'Time' coordinate of xarray data sets must be consistent 
with the MPAS calendar <br>
Date last modified: 2017/02/09 <br>
Contributors: Xylar Asay-Davis
</h2>

For all data sets used in the analysis, the 'Time' coordinate must represent dates 
on the appropriate MPAS calendar, either 'gregorian' or 'gregorian_noleap', depending
on the namelist option 'config_calendar_type'.  There must be ways of mathematically
manipulating times (e.g. adding/subtracting offsets and figuring out the amount of time
between two dates) and of making plots that are consistent with these calendars.

<h2> Requirement: The 'Time' coordinate of xarray data sets must support at least years 
0001 and 9999, and preferably any conceivable value<br>
Date last modified: 2017/02/16 <br>
Contributors: Xylar Asay-Davis
</h2>

For all data sets used in the analysis, the 'Time' coordinate must, at a minimum, 
support years between 0001 and 9999 (the range of `datetime.datetime`) and preferably
a broader range.


<h1> Algorithmic Formulations (optional) </h1>

<h2> Design solution: The 'Time' coordinate of xarray data sets must be consistent
with the MPAS calendar <br>
Date last modified: 2017/02/11 <br>
Contributors: Xylar Asay-Davis, Phillip J. Wolfram
</h2>

The proposed solution represents time in `xarray.DataSet` objects as the number of
days since the reference date 0001-01-01.  
This is reasonable because the smallest unit of time output in MPAS components is
seconds (and unlikely to ever be shorter than ms).  We note that a date specified
as a 64-bit float has a precision high enough to represent seconds for dates up
to +/- 100 million years: 
```python
>>> import sys
>>> 1./(sys.float_info.epsilon*365*24*60*60)
142808207.36207813
```
We should have no trouble representing any number we might want (including paleo
timescales) with this system.

For purposes of performing mathematical operations and plotting dates, these
values will be converted to `datetime.datetime` objects (via the proposed
`days_to_datetime` utility function) and back (via the proposed 
`datetime_to_days`).

The conversion operations within `datetime_to_days` and `days_to_datetime` will be
performed with the calendar-aware functions `netCDF4.date2num` and 
`netCDF4.num2date`, respectively.  Both functions will support lists/arrays of dates
(for efficiency and simplicity of calling code) in addition to single values.

Curve ploting can be supported with `matplotlib.pyplot.plot_date`, which takes a date
of exactly the format used here (days since 0001-01-01).  The compatibility with `plot_date`
was part of the reason for choosing this format for the date.

<h2> Design solution: The 'Time' coordinate of xarray data sets must support at least years 
0001 and 9999, and preferably any conceivable value<br>
Date last modified: 2017/02/09 <br>
Contributors: Xylar Asay-Davis
</h2>

Same as above. In theory, the use of days since 0001-01-01 would allow any year
to be supported, not just the range from 0001 to 9999.  However, the conversions
to `datetime.datetime` objects for mathematical manipulation will constrain
the dates to be between `datetime.min` (0001-01-01) and `datetime.max` (9999-12-31).


<h1> Design and Implementation </h1>

<h2> Implementation: The 'Time' coordinate of xarray data sets must be consistent
with the MPAS calendar <br>
Date last modified: 2017/02/16 <br>
Contributors: Xylar Asay-Davis
</h2>

The proposed implementation is on the branch 
[xylar/generalize_calendar](https://github.com/xylar/MPAS-Analysis/tree/generalize_calendar)

A helper funciton, `mpas_xarray._parse_dataset_time`, computes times as days since 
0001-01-01, and serves as a replacement for `mpas_xarray._get_datetimes`. 

**Note: the current implementation breaks the convention that `mpas_xarray` remains
separate from the rest of MPAS-Analyis by using 3 functions from `timekeeping.utility`
in `mpas_xarray`:**
```python
from ..timekeeping.utility import string_to_days_since_date, \
    days_to_datetime, datetime_to_days
```
**This violates the first requirement in the 
[Design Document: Moving variable mapping out of mpas_xarray](https://github.com/xylar/MPAS-Analysis/blob/design_doc_variable_mapping_reorg/design_docs/variable_mapping_reorg.md).
I am open to alternative solutions for keeping `mpas_xarray` separate from the rest
of analysis but these 3 functions do not conceptually belong in `mpas_xarray`.  The
problem is exacerbated by the fact that there are analysis-specific functions in
`timekeeping`, meaning that this cannot easily be made a submodule of `mpas_xarray`
(nor would this make very much logical sense).  Having 2 `timekeeping` modules, one
for `mpas_xarray` and one for MPAS-Analysis, seems unnecessarily confunsing.**

The functions `generalized_reader.open_multifile_dataset` and
`mpas_xarray.open_multifile_dataset` have been updated to use this method for parsing
times.  This involves removing the `year_offset` argument and adding an optional 
`simulation_start_time` argument for supplying a date to use to convert variables
like `daysSinceStartOfSim` to days since 0001-01-01.

An example of opening a data set and manipulating times withe the new approach in
the OHC script is:
```python
from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime, string_to_datetime
...
def ohc_timeseries(config, streamMap=None, variableMap=None):
...
    simulationStartTime = get_simulation_start_time(streams)
...
    ds = open_multifile_dataset(file_names=file_names,
                                calendar=calendar,
                                simulation_start_time=simulation_start_time,
                                time_variable_name='Time',
                                variable_list=variable_list,
                                variable_map=variableMap,
                                start_date=startDate,
                                end_date=endDate)

    timeStart = string_to_datetime(startDate)
    timeEnd = string_to_datetime(endDate)

    # Select year-1 data and average it (for later computing anomalies)
    timeStartFirstYear = string_to_datetime(simulation_start_time)
    if timeStartFirstYear < timeStart:
        startDateFirstYear = simulation_start_time
        firstYear = int(startDateFirstYear[0:4])
        endDateFirstYear = '{:04d}-12-31_23:59:59'.format(firstYear)
        filesFirstYear = streams.readpath(streamName,
                                          startDate=startDateFirstYear,
                                          endDate=endDateFirstYear,
                                          calendar=calendar)
        dsFirstYear = open_multifile_dataset(
            file_names=filesFirstYear,
            calendar=calendar,
            simulation_start_time=simulation_start_time,
            time_variable_name='Time',
            variable_list=variable_list,
            variable_map=variableMap,
            start_date=startDateFirstYear,
            end_date=endDateFirstYear)
    else:
        dsFirstYear = ds
        firstYear = timeStart.year

    timeStartFirstYear = date_to_days(year=firstYear, month=1, day=1,
                                      calendar=calendar)
    timeEndFirstYear = date_to_days(year=firstYear, month=12, day=31,
                                    hour=23, minute=59, second=59,
                                    calendar=calendar)

    dsFirstYear = dsFirstYear.sel(Time=slice(timeStartFirstYear,
                                             timeEndFirstYear))

    meanFirstYear = dsFirstYear.mean('Time')
...
    yearStart = days_to_datetime(ds.Time.min()).year
    yearEnd = days_to_datetime(ds.Time.max()).year
    timeStart = date_to_days(year=yearStart, month=1, day=1,
                             calendar=calendar)
    timeEnd = date_to_days(year=yearEnd, month=12, day=31,
                           calendar=calendar)
                           
    if preprocessedReferenceRunName != 'None':
        print '  Load in OHC from preprocessed reference run...'
        inFilesPreprocessed = '{}/OHC.{}.year*.nc'.format(
            preprocessedInputDirectory, preprocessedReferenceRunName)
        dsPreprocessed = open_multifile_dataset(
            file_names=inFilesPreprocessed,
            calendar=calendar,
            simulation_start_time=simulation_start_time,
            time_variable_name='xtime')
        yearEndPreprocessed = days_to_datetime(dsPreprocessed.Time.max()).year
...
```

The `replicate_cycles` function in `sea_ice.timeseries` has been a particular
challenge with the existing calendar.  Here is that function with the new 'Time'
coordinate:
```python
def replicate_cycle(ds, dsToReplicate, calendar):
    dsStartTime = days_to_datetime(ds.Time.min(), calendar=calendar)
    dsEndTime = days_to_datetime(ds.Time.max(), calendar=calendar)
    repStartTime = days_to_datetime(dsToReplicate.Time.min(),
                                    calendar=calendar)
    repEndTime = days_to_datetime(dsToReplicate.Time.max(),
                                  calendar=calendar)

    repSecondTime = days_to_datetime(dsToReplicate.Time.isel(Time=1),
                                     calendar=calendar)

    period = (MpasRelativeDelta(repEndTime, repStartTime) +
              MpasRelativeDelta(repSecondTime, repStartTime))

    startIndex = 0
    while(dsStartTime > repStartTime + (startIndex+1)*period):
        startIndex += 1

    endIndex = 0
    while(dsEndTime > repEndTime + (endIndex+1)*period):
        endIndex += 1

    dsShift = dsToReplicate.copy()

    times = days_to_datetime(dsShift.Time, calendar=calendar)
    dsShift.coords['Time'] = ('Time',
                              datetime_to_days(times + startIndex*period,
                                               calendar=calendar))
    # replicate cycle:
    for cycleIndex in range(startIndex, endIndex):
        dsNew = dsToReplicate.copy()
        dsNew.coords['Time'] = ('Time',
                                datetime_to_days(times + (cycleIndex+1)*period,
                                                 calendar=calendar))
        dsShift = xr.concat([dsShift, dsNew], dim='Time')

    return dsShift
```

<h2> Implementation: The 'Time' coordinate of xarray data sets must support at least years 
0001 and 9999, and preferably any conceivable value<br>
Date last modified: 2017/02/09 <br>
Contributors: Xylar Asay-Davis
</h2>

Same as above.

<h1> Testing </h1>

<h2> Testing and Validation: The 'Time' coordinate of xarray data sets must be consistent
with the MPAS calendar <br>
Date last modified: 2017/02/11 <br>
Contributors: Xylar Asay-Davis
</h2>
In [xylar/generalize_calendar](https://github.com/xylar/MPAS-Analysis/tree/generalize_calendar), 
unit testing has been added for `timekeeping` and `mpas_xarray` that checks both the `gregorian` 
and `gregorian_noleap` calendars under simple test conditions.  However, we have no data sets 
that test `gregorian`, so we have a somewhat limited ability to test this calendar option.
Fortunately, there are also no immediate plans to run with `gregorian`.

I will make sure all tests with config files in the `configs/lanl` and `configs/edison` 
directories produce bit-for-bit results with the current `develop`.

<h2> Testing and Validation: The 'Time' coordinate of xarray data sets must support at least years 
0001 and 9999, and preferably any conceivable value<br>
Date last modified: 2017/02/11 <br>
Contributors: Xylar Asay-Davis
</h2>

Unit tests have been added to ensure that dates both close to 0001-01-01 and typical
calendar dates (e.g. 2017-01-01) function as expected.

@akturner's MPAS-SeaIce run with real dates (mentioned in 
[#81](https://github.com/MPAS-Dev/MPAS-Analysis/issues/81))  has been successfully
run with the proposed approach.  This run started in 1958, and had presented a problem
for MPAS-Analysis with the previous calendar.
