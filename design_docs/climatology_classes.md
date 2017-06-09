<h1> Title: Changing climatology module to use classes <br>
Xylar Asay-Davis <br>
date: 2017/06/09 <br>
</h1>
<h2> Summary </h2>
After recent changes to MPAS-Analysis, tasks are represented by classes
(meaning that there are objects for each task that carry their state with
 them).  In making modifications to these tasks, it became clear to me that
there was a lot of redundant code in individual tasks related to computing
climatologies, and that this redundant code could be moved into the
climatology module itself if climatologies, too, became classes.  This would
mean that climatology objects could carry around their state (e.g. a
climatology data set, a remapper, a remapped data set, the task for which the
climatology is being computed). The move to climatology objects would mean that
a task can construct a climatology object, passing in the necessary information
to automatically create descriptors of the source and destination (comparison)
grids, choose appropriate file names, create a remapper (if necessary) and/or
load a remapped data set (if one has already been cached) without individual
tasks needing to handle these details.  This design will be considered
successful if the code needed to compute and optionally remap a climatology
is significantly reduced.


<h1> Design and Implementation </h1>

I am skipping the Requirements section because this is largely a change in
implementation that does not alter the requirements of the earlier climatology
design.

The design has been implemented in the branch
[xylar/climatology_caching_skip_opening_ds](https://github.com/xylar/MPAS-Analysis/tree/climatology_caching_skip_opening_ds)
and should be merged through pull request
[#204](https://github.com/MPAS-Dev/MPAS-Analysis/pull/204)


<h2> Implementation: `Climatology` base class <br>
Date last modified: 2017/06/09 <br>
Contributors: Xylar Asay-Davis
</h2>

The `Climatology` base class can be used directly to compute a climatology from
a given `xarray` data set.  A seasonal climatology can be computed via `compute`, 
either supplying the name of a season via `monthNames` or a list of months via 
the `monthValues` argument.  A set of monthly climatologies (1 for each month) can
be computed via `compute_monthly`.  `Climatology` also implements 2 methods, 
`create_remapper` and `remap_and_write`, that assist in remapping the climatology
data set onto a comparison grid.  These require additional attributes to be defined
(see below) so are typically not used directly by `Climatology` objects but rather
by objects of one of the derived classes (`MpasClimatology` or 
`ObservationsClimatology`).

The class is structured as follows:
```python
class Climatology(object):
    """
    A class for computing climatologies from a monthly mean data set.
    """

    def __init__(self, task, monthNames=None, comparisonGrid=None):
        """
        Create a new climatology.
        """

    def create_remapper(self, mappingFileSection, mappingFileOption,
                        mappingFilePrefix, method):
        """
        Creates an attribute ``remapper``, a ``Remapper`` object, that can be
        used to remap from source files or data sets to corresponding data sets
        on the comparison grid.
        """

    def compute(self, ds, monthValues=None, maskVaries=True):
        """
        Compute a monthly, seasonal or annual climatology data set from a data
        set.  The mean is weighted but the number of days in each month of
        the data set, ignoring values masked out with NaNs.  If the month
        coordinate is not present, a data array ``month`` will be added based
        on ``Time`` and the provided calendar.
        """

    def compute_monthly(self, ds, maskVaries=True):  # {{{
        """
        Compute monthly climatologies from a data set.  The mean is
        weighted by the number of days in each month of the data set,
        ignoring values masked out with NaNs.  If the month coordinate is
        not present, a data array ``month`` will be added based on
        ``Time``.
        """
        
    def remap_and_write(self, useNcremap=None):
        """
        Remap the climatology data set produced by a call to ``compute`` or
        ``cache``, write the result to an output file, and return the remapped
        data set.
        """
```

The class contains a number of attributes:
- `task` : `AnalysisTask` object <br>
        An analysis task for which the climatology is needed. `task` is
        used to get config options, define the component name, map streams,
        namelists and variables, etc.

- `config` :  instance of `MpasAnalysisConfigParser` <br>
        Contains configuration options

- `calendar` : {'gregorian', 'gregorian_noleap'} <br>
        the name of the calendar

- `monthNames` : str  <br>
        The months that make up the climatology, used in constructing the
        names of cache files. If provided, `monthNames` should be one of
        the keys of `monthDictionary` in the `constants` module.

- `monthValues` : list of int <br>
        A list of integer months that make up the season in `monthNames`,
        taken from `monthDictionary` in the `constants` module.  The
        entries are sorted in ascending order.

- dataSet : `xarray.Dataset` or `xarray.DataArray` object <br>
        A data set containing the climatology (before remapping, if
        applicable).  Available only after calling the `compute` method.

Several more attributes are only defined if remapping is to be performed
(indicted by the value of `comparisonGrid` passed to `__init__`):
- `comparisonDescriptor` : `MeshDescriptor` object <br>
        If the name of a  comparison grid name was supplied, this is an object
        describing that comparison grid (e.g. for remapping).  If no comparison
        grid name was supplied, this is `None`.
- `remapper` : `Remapper` object <br>
        A remapper between the source and comparison grids.  Available only
        after calling the `create_remapper` method.
- `remappedDataSet` : `xarray.Dataset` or `xarray.DataArray` object <br>
        A data set containing the remapped climatology.  Available only after
        calling the `remap_and_write` method.

Several more attributes must be defined by a child class if remapping is to 
be performed:

- `sourceDescriptor` : `MeshDescriptor` object <br>
        A descriptor of the source grid or mesh used for remapping.  This
        attribute should be set by one of the subclasses of `Climatology` if
        remapping will be performed.

- climatologyFileName : str <br>
        The name of the file where `dataSet` (the climatology data set before
        remapping) should be stored.  This attribute should be set by one of
        the subclasses of `Climatology` if remapping will be performed.

- remappedFileName : str <br>
        The name of the file where `remappedDataSet` (the climatology data
        set after remapping) should be stored.  This attribute should be set by
        one of the subclasses of `Climatology` if remapping will be
        performed.

Example usage of `Climatology` to compute a seasonal climatology for
January, February and March from within an analysis task using a data
set `ds` would look something like this:
```python
from ..shared.climatology import Climatology
...
climatology = Climatology(task=self, monthNames='JFM')
climatology.compute(ds=ds)
# do computations or plotting on climatology.dataSet
```

Computing a monthly climatology from `ds` would look like this:
```python
from ..shared.climatology import Climatology
...
climatology = Climatology(task=self)
climatology.compute_monthly(ds=ds)
# do computations or plotting on climatology.dataSet
```

<h2> Implementation: `MpasClimatology` class <br>
Date last modified: 2017/06/09 <br>
Contributors: Xylar Asay-Davis
</h2>

The `MpasClimatology` class implements functionality for handling data sets on
MPAS meshes.  In the process of constructing a climatology object, file names
for cache files are automatically determined from the field being cached, the
name of the MPAS mesh, the years in the climatology, etc. If the name of a
comparison grid is provided (currently `latlon` is the only supported option,
corresponding to the default global lat-lon grid), the climatology object will
handle remapping data or loading data from an existing remapped cache file.

The class is structured as follows:
```python
class MpasClimatology(Climatology):  # {{{
    """
    A class for computing climatologies from an MPAS monthly mean data set.
    """
    def __init__(self, task, fieldName, monthNames, streamName,
                 meshFileName=None, comparisonGrid=None,
                 mappingFileSection=None, mappingFileOption=None,
                 mappingFilePrefix=None, method=None):  # {{{
        """
        Create a new object for computing climatologies from an MPAS monthly
        mean data set.

        If a ``comparisonGrid`` is provided, the climatology object also
        either loads data from a cache file containing remapped data (if
        one already exists) or creates a remapper that can be used later
        on (via the ``remap_and_write`` method) to compute the remapped
        climatology.
        """

def cache(self, openDataSetFunc, printProgress=False):  # {{{
        '''
        Cache NetCDF files for each year of an annual climatology, and then use
        the cached files to compute a climatology for the full range of years.
        The start and end years of the climatology are taken from ``config``,
        and are updated in ``config`` if the data set ``ds`` doesn't contain
        this full range.

        Note: only works with climatologies where the mask (locations of
        ``NaN`` values) doesn't vary with time.
        """
```

Other than setting up file names and the ``remapper`` needed for caching and
remapping climatologies, the primary functionality added by `MpasClimatology`
is the `cache` method.  This method acts somewhat like `compute` but it takes
a function `openDataSetFunc` as an argument instead of an existing data set.
This function is used to open part of a data set during the caching process.
The function must have arguments `inputFileNames`, `startDate` and `endDate`,
and might look something like the following, implemented as a private
method of a task:
```python
    def _open_mpas_dataset_part(self, inputFileNames, startDate,
                                endDate):  # {{{
        """
        Open part of a data set between the given start and end date, used
        to cache a climatology of the data set.

        Parameters
        ----------
        inputFileNames : list of str
            File names in the multifile data set to open

        startDate, endDate : float
            start and end date to which to crop the Time dimension (given in
            days since 0001-01-01)
        """
        ds = open_multifile_dataset(
            fileNames=inputFileNames,
            calendar=self.calendar,
            config=self.config,
            simulationStartTime=self.simulationStartTime,
            timeVariableName='Time',
            variableList=['sst'],
            variableMap=self.variableMap,
            startDate=startDate,
            endDate=endDate)
        return ds  # }}}
```
Note that, before the data set is returned, one could perform further manipulation
such as adding additional data arrays or coordinates.

`MpasClimatology` also inherits all the methods of `Climatology` (`create_remapper`,
`compute`, `compute_monthly` and `remap_and_write`) as well as all of its attributes.
Note, however, that `creat_remapper` gets called automatically when an `MpasClimatology`
is created with a comparison grid name given as an argument, so this method does not
need to be called directly.

During construction, 3 attributes needed for performing remapping are defined:

- `sourceDescriptor` : `MpasMeshDescriptor` object <br>
        A descriptor of the source MPAS mesh used for remapping (defined only if
        a comparison grid name is provided)

- `climatologyFileName` : str <br>
        The name of the cache file where `dataSet` (the climatology data set before
        remapping) should be stored.

- `remappedFileName` : str <br>
        The name of the cache file where `remappedDataSet` (the climatology data
        set after remapping) should be stored (defined only if a comparison grid 
        name is provided)

There are also several attributes that are specific to `MpasClimatology` (not
inherited from or used by `Climatology`):

- `fieldName` : str <br>
        The name of the field for which the climatology is being computed,
        used in constructing the names of cache files.

- `streamName` : str <br>
        The name of the stream from which the climatology data set will
        be read, used to cache times and update their bounds in the
        configuration parser.

- `startYear`, `endYear` : int <br>
        The start and end years of the climatology

- `climatologyDirectory` : str <br>
        The directory where climatologies will be cached

- `climatologyPrefix` : str <br>
        The prefix (including full path) on climatology cache files.  A call
        to the `cache` method may produce several intermediate cache files
        as well as a "final" cache file with the seasonal average between
        `startYear` and `endYear`.

Example usage of `MpasClimatology` to compute an annual climatology of SST 
(sea-surface temperature) from within an analysis task would look something
like this:
```python
from ..shared.climatology import MpasClimatology
...
restartFileName = self.runStreams.readpath('restart')[0]
climatology = MpasClimatology(
        task=self,
        fieldName='sst',
        monthNames='ANN',
        streamName='timeSeriesStats',
        meshFileName=restartFileName,
        comparisonGrid='latlon',
        mappingFileSection='climatology',
        mappingFileOption='mpasMappingFile',
        mappingFilePrefix='map',
        method='bilinear')
# if the remapped data set has already been cached, it will
# have been read into climatology.remappedDataSet
if climatology.remappedDataSet is None:
    # cache the climatology (computing incremental climatologies along
    # the way)
    climatology.cache(
            openDataSetFunc=self._open_mpas_dataset_part,
            printProgress=True)
    # Now that we have a climatology (stored in climatology.dataSet), we
    # can remap it and write it out to the cache file
    climatology.remap_and_write()
# Either way, we now should have a valid remapped climatology in 
# climatology.remappedDataSet that we can perform computations on and plot.
```

A simpler example that creates an annual climatology but which does not perform
remapping is as follows:
```python
from ..shared.climatology import MpasClimatology
...
climatology = MpasClimatology(
        task=self,
        fieldName='sst',
        monthNames='ANN',
        streamName='timeSeriesStats')
climatology.cache(
        openDataSetFunc=self._open_mpas_dataset_part,
        printProgress=True)
# We now should have a valid annual climatology in climatology.dataSet that we can
# perform computations on and plot.
```

<h2> Implementation: `ObservationsClimatology` class <br>
Date last modified: 2017/06/09 <br>
Contributors: Xylar Asay-Davis
</h2>

The `ObservationsClimatology` class implements functionality for handling 
observational data sets.  In the process of constructing a climatology object, 
file names for cache files are automatically determined from the field being 
cached, the name of the grid, etc. If the name of a comparison grid is provided
(currently `latlon` is the only supported option, corresponding to the default 
global lat-lon grid), the climatology object will handle remapping data or 
loading data from an existing remapped cache file.

The class is structured as follows:
```python
class ObservationsClimatology(Climatology):  # {{{
    """
    A class for computing climatologies from an obsevational data set.
    """
    def __init__(self, task, fieldName, monthNames,
                 obsGridDescriptor, comparisonGrid=None,
                 mappingFileSection=None, mappingFileOption=None,
                 mappingFilePrefix=None, method=None):  # {{{
        """
        Create a new object for creating climatologies from observational data
        sets.
        """
```

A `LatLonGridDescriptor` or `ProjectionGridDescriptor` object must be supplied as 
an input argument (`obsGridDescriptor`) when constructing the climatology object. 
This descriptor can either be created (e.g. via the `LatLonGridDescriptor.create` method) 
by supplying arrays describing the coordinates or read from a file or opened data set(e.g. 
via the `LatLonGridDescriptor.read` method).  Here is an example where we manipulate the
data set before using it to create a `LatLonGridDescriptor`, then compute and remap an
annual climatology of mixed-layer depth (MLD):
```python
from ..shared.grid import LatLonGridDescriptor
...
dsObs = xr.open_mfdataset('holtetalley_mld_climatology.nc')
dsObs.rename({'lat': 'latCoord', 'lon': 'lonCoord'}, inplace=True)
dsObs.rename({'iLAT': 'lat', 'iLON': 'lon'}, inplace=True)
# set the coordinates now that the dimensions have the same names
dsObs.coords['lat'] = dsObs['latCoord']
dsObs.coords['lon'] = dsObs['lonCoord']

obsDescriptor = LatLonGridDescriptor.read(ds=dsObs)
climatology = ObservationClimatology(
        task=self,
        fieldName='mld',
        monthNames='ANN',
        obsGridDescriptor=obsDescriptor,
        comparisonGrid='latlon',
        mappingFileSection='oceanObservations',
        mappingFileOption='mldClimatologyMappingFile',
        mappingFilePrefix='map_obs_mld',
        method='bilinear')
if climatology.remappedDataSet is None:
    climatology.compute(ds=dsObs)
    climatology.remap_and_write()
# Either way, we now should have a valid remapped climatology in 
# climatology.remappedDataSet that we can perform computations on and plot.
```

`ObservationsClimatology` does not support sophisticated caching like
`MpasClimatology` because it is assumed that monthly climatologies of
observations have been precomputed.

`ObservationsClimatology` inherits all the methods of `Climatology` (`create_remapper`,
`compute`, `compute_monthly` and `remap_and_write`) as well as all of its attributes.
Note that `creat_remapper` gets called automatically when an `ObservationsClimatology`
is created with a comparison grid name given as an argument, so this method does not
need to be called directly.

During construction, 3 attributes needed for performing remapping are defined:

- `sourceDescriptor` : `MeshDescriptor` object <br>
        A descriptor of the source grid used for remapping (defined only if
        a comparison grid name is provided)

- `climatologyFileName` : str <br>
        The name of the cache file where `dataSet` (the climatology data set before
        remapping) should be stored.

- `remappedFileName` : str <br>
        The name of the cache file where `remappedDataSet` (the climatology data
        set after remapping) should be stored (defined only if a comparison grid 
        name is provided)

There are also 2 attributes that are specific to `ObservationsClimatology` (not
inherited from or used by `Climatology`):

- `fieldName` : str <br>
        The name of the field for which the climatology is being computed,
        used in constructing the names of cache files.

- `climatologyDirectory` : str <br>
        The directory where climatologies will be cached


<h1> Testing and Validation</h1>
<h2> Date last modified: 2017/06/09 <br>
Contributors: Xylar Asay-Davis
</h2>

All plots will be tested to ensure they are bit-for-bit identical to
those produced by `develop` for all tests defined in the `configs/edison`
and `configs/lanl` directories.
