<h1> Title: Caching Times from Multifile Data Sets <br>
Xylar Asay-Davis <br>
date: 2017/06/10 <br>
</h1>
<h2> Summary </h2>
Currently, we cache climatologies and time series computed from mutlifile data
sets.  However, each time we perform the analysis, we open all the files in the
data set (once per analysis task) to determine what times each file contains
and whether this data has already been cached.  A mechanism is needed to cache
the times in each multifile data set that avoids opening files multiple times
if climatology or time series data has already been computed and cached. This
allows the user to change the bounds of a climatology or time series and to
rerun the analysis without having to recompute data that has already been
cached. It also means that the analysis will automatically only access new
files if a simulation is underway and new data becomes available that should
be included in a time series or climatology.

Success of this design will mean that MPAS data files in a multifile data set
are not accessed only the first time MPAS-Analysis is run on those files
(unless the files have been modified).

<h1> Requirements </h1>

<h2> Requirement: Cache times from multifile data sets <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

A mechanism should exist for storing the times in a multifile data set so that
the files in that data set do not need to be opened just to determine which
times they contain.  Since multiple tasks use the same multifile data set in
many cases, the times in that data set should be cached only once and then used
for multiple tasks.

<h2> Requirement: Use caching times when updating climatologies and time 
series <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

A mechanism should exist for reading the cached times and using them to inform
the code that updates cached climatologies and time series.  

<h1> Algorithmic Formulations </h1>
<h2> Design solution: Cache times from multifile data sets <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

The cached times will be stored in a dictionary of dictionaries.  The outer
dictionary will have keys that are the file names in the data set.  Each inner
dictionary contains:
- `times`: the times (in days since 0001-01-01) contained in the file
- `years`, `months`: the year and month corresponding to each entry in `times`
- `daysInMonth`: the number of days over which the data has been averaged in
   each entry in the data file (the days in a month, since this is assumed to
   be monthly averaged data)
- `dateModified`: the date the file was last modified (used to update the 
   cached times if needed)
This dictionary will be stored in a cache file (in python "pickle" format,
since NetCDF does not have a good way of storing python dictionaries).

To be compatible with parallel tasks, caching will be performed by an analysis
task (defined in class `CacheDatasetTimesTask`) that will be a prerequisite
of any tasks using that multifile data set.  Details of 
[prerequisite tasks](prerequisite_tasks.md) are discussed in a separate design
document.   A task defined by `CacheDatasetTimesTask` will take a component
name, a stream name and a list of config sections containing entries 
`startYear` and `endYear` (currently `climatology`, `timeSeries` and `index`)
that are used to determine which files should be included in the time caching
(i.e. no need to cache times from files that we don't currently want to plot).

<h1> Algorithmic Formulations </h1>
<h2> Design solution: Use caching times when updating climatologies and time 
series <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

As part of caching climatologies, the dictionary of cached times will be read
in from a cache file.  Then, the code will determine if any files contain times
that are not already included in cached climatologies and compute (or update)
the cached files accordingly.  Intermediate climatologies (e.g. yearly 
climatologies if `yearsPerCacheFile = 1` ) will be recomputed if they are
"incomplete" (e.g. if a yearly climatology doesn't contain all 12 months).  Any
cache files that are "complete" will not be recomputed, meaning that the
corresponding files in the multifile data set should not be accessed at all.

<h1> Design and Implementation </h1>

The design has been implemented in the branch
[xylar/climatology_caching_skip_opening_ds](https://github.com/xylar/MPAS-Analysis/tree/climatology_caching_skip_opening_ds)
and is expected to be merged through pull request
[#204](https://github.com/MPAS-Dev/MPAS-Analysis/pull/204)

<h2> Implementation: Cache times from multifile data sets <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

Climatologies are cached through the `cache_multifile_dataset_times` method
of `AnalysisTask`:
```python
def cache_multifile_dataset_times(self, inFileNames, streamName,
                                  timeVariableName='Time'):
    """
    Creates a cache file of the times in each file of a multifile data set.
    This is useful when caching climatologies and time series as a
    simulation evolves, since files that have already been processed will
    not need to be opened to find out which times they contain.

    Parameters
    ----------
    inFileNames : list of str
        A list of file paths to read

    streamName : str
        The name of a stream, used to build the name of the cache file

    timeVariableName : string, optional
        The name of the time variable (typically 'Time' if using a
        variableMap or 'xtime' if not using a variableMap)
    """
```
The stream name is typically `timeSeriesStats` but could be any stream
(possibly renamed through the streams map to be compatible with multiple
MPAS versions). The input files can be determined by calling the 
`get_input_file_names` method of `AnalysisTask` with the name of the stream 
and either a start and end date or a config section from which to read the 
start and end date.

```python
def get_input_file_names(self, streamName,
                         startDate=None, endDate=None,
                         startAndEndDateSection=None):
    '''
    Get a list of input files corresponding to the given stream and
    optionally bounded by the start and end dates found in the given
    section of the config file.

    Parameters
    ----------
    streamName : str
        The name of a stream to check.  If ``self.streamMap`` is defined,
        the streamName will be mapped to the corresponding name in the
        streams file

    startDate, endDate : float, optional
        start and end date to use in determining which files to include in
        the list

    startAndEndDateSection : str, optional
        If ``startDate`` and ``endDate`` arguments are not supplied, the
        name of a section in the config file containing ``startDate`` and
        ``endDate`` options to use instead. ``startAndEndDateSection`` is
        typically one of ``climatology``, ``timeSeries`` or ``index``.
    '''
```

`CacheDatasetTimesTask`, the analysis task for caching times, is in the file
`mpas_analysis/shared/cache_dataset_times_task.py`.  The constructor is as
follows:
```python
class CacheDatasetTimesTask(AnalysisTask):
    '''
    A task for caching the times in a multifile data sets of an MPAS analysis
    member for later use.  Since analysis member may be used by multiple tasks
    (indeed the ``timeSeriesStats`` is currently used by *all* tasks), it is
    important that this time information gets processed once before all other
    tasks get run in parallel.
say-Davis
    '''

    def __init__(self, config, componentName, streamName,
                 startAndEndDateSections, namelistOption=None):
        '''
        Construct an analysis task for caching the times in the multifile
        data set in the given component and stream.  The name of the task
        includes the component and stream name.  For example, if
        ``component='ocean'`` and ``streamName='timeSeriesStats``, then the
        task name is ``cacheOceanTimeSeriesStatsTimes``.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        componentName :  {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        streamName : str
            The name of the stream from which the climatology data set will
            be read, used to cache times and update their bounds in the
            configuration parser.

        startAndEndDateSections : list, {'climatology', 'timeSeries', 'index'}
            The name of sections in the config file containing ``startDate``
            and ``endDate`` options as a list.

        namelistOption : str, optional
            The name of a namelist option (e.g.
            ``config_am_timeseriesstatsmonthly_enable``) that should be set to
            true of the required analysis member has been enabled.  If this
            option is ``None`` (the default), no check is performed

        '''
```
Note that the analysis task is **not** named `cacheDatasetTimesTask`, as would
typically be the case base don the name of the class.  Instead, because this
class can be used to create multiple analysis tasks, each with a different
component and stream name, the name of the task is 
`cache<Component><StreamName>Times` (e.g. `cacheOceanTimeSeriesStatsTimes`),
which should lead to time-caching analysis task having a unique name.

<h2> Implementation: Use caching times when updating climatologies and time 
series <br>
Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>

The cached times are read back using the same method of `AnalysisTask`, 
`cache_multifile_dataset_times`, that was used to compute them and read them
in the first place. If no files have been modified since the last caching,
the cache file containing the times is simply read in and the dictionary is
returned, as in this example from the `MpasClimatology` class:
```python
...
startDate = self.config.get('climatology', 'startDate')
endDate = self.config.get('climatology', 'endDate')
startDate = string_to_days_since_date(dateString=startDate,
                                      calendar=self.calendar)
endDate = string_to_days_since_date(dateString=endDate,
                                    calendar=self.calendar)

inFileNames = self.task.get_input_file_names(
    self.streamName, startAndEndDateSection='climatology')

fullTimeCache = self.task.cache_multifile_dataset_times(
    inFileNames, self.streamName, timeVariableName='Time')

# find only those cached times between starDate and endDate
timeCache = OrderedDict()
times = []
for fileName in fullTimeCache:
    localTimes = fullTimeCache[fileName]['times']
    mask = numpy.logical_and(localTimes >= startDate,
                             localTimes < endDate)
    if numpy.count_nonzero(mask) == 0:
        continue

    times.extend(list(localTimes[mask]))
    timeCache[fileName] = fullTimeCache[fileName]

...
```
timeCache is then used to determine which files are needed to make up each
intermediate climatology during climatology caching an only those files for
for which cached climatologies are "incomplete" are accessed.


<h1> Testing and Validation </h1>
<h2> Date last modified: 2017/06/10 <br>
Contributors: Xylar Asay-Davis
</h2>
All plots will be tested to ensure they are bit-for-bit identical to
those produced by `develop` for all tests defined in the `configs/edison`
and `configs/lanl` directories.

In addition, a test will be run where with climatologies computed over years
20-50 or a long data set, then over 20-60, making sure that only files covering
50-60 are accessed the second time
