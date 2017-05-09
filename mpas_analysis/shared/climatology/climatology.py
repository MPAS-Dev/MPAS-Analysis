"""
Functions for creating climatologies from monthly time series data

Authors
-------
Xylar Asay-Davis
"""

import xarray as xr
import os
import numpy
import warnings
from collections import OrderedDict
import numbers

from ..constants import constants

from ..timekeeping.utility import string_to_days_since_date, \
    add_years_months_days_in_month

from ..io import build_config_full_path, make_directories
from ..io.utility import fingerprint_generator

from ..interpolation import Remapper
from ..grid import MpasMeshDescriptor, LatLonGridDescriptor, \
    ProjectionGridDescriptor


def get_lat_lon_comparison_descriptor(config):  # {{{
    """
    Get a descriptor of the lat/lon comparison grid, used for remapping and
    determining the grid name

    Parameters
    ----------
    config :  instance of ``MpasAnalysisConfigParser``
        Contains configuration options

    Returns
    -------
    descriptor : ``LatLonGridDescriptor`` object
        A descriptor of the lat/lon grid

    Authors
    -------
    Xylar Asay-Davis
    """
    climSection = 'climatology'

    comparisonLatRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLatitude)
    comparisonLonRes = config.getWithDefault(climSection,
                                             'comparisonLatResolution',
                                             constants.dLongitude)

    nLat = int((constants.latmax-constants.latmin)/comparisonLatRes)+1
    nLon = int((constants.lonmax-constants.lonmin)/comparisonLonRes)+1
    lat = numpy.linspace(constants.latmin, constants.latmax, nLat)
    lon = numpy.linspace(constants.lonmin, constants.lonmax, nLon)

    descriptor = LatLonGridDescriptor.create(lat, lon, units='degrees')

    return descriptor  # }}}


class Climatology(object):  # {{{
    """
    A class for computing climatologies from a monthly mean data set.

    Attributes
    ----------
    task : ``AnalysisTask`` object
        An analysis task for which the climatology is needed. ``task`` is
        used to get config options, define the component name, map streams,
        namelists and variables, etc.

    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    calendar : {'gregorian', 'gregorian_noleap'}
        the name of the calendar

    monthNames : str
        The months that make up the climatology, used in constructing the
        names of cache files. If provided, ``monthNames`` should be one of
        the keys of ``monthDictionary`` in the ``constants`` module.

    monthValues : list of int
        A list of integer months that make up the season in ``monthNames``,
        taken from ``monthDictionary`` in the ``constants`` module.  The
        entries are sorted in ascending order.

    sourceDescriptor : ``MeshDescriptor`` object
        A descriptor of the source grid or mesh used for remapping.  This
        attribute should be set by one of the subclasses of ``Climatology`` if
        remapping will be performed.

    comparisonDescriptor : ``MeshDescriptor`` object
        If the name of a  comparison grid name was supplied, this is an object
        describing that comparison grid (e.g. for remapping).  If no comparison
        grid name was supplied, this is ``None``.

    remapper : ``Remapper`` object
        A remapper between the source and comparison grids.  Available only
        after calling the ``create_remapper`` method.

    dataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the climatology (before remapping, if
        applicable).  Available only after calling the ``compute`` or
        ``compute_monthly`` method.

    climatologyFileName : str
        The name of the file where ``dataSet`` (the climatology data set before
        remapping) should be stored.  This attribute should be set by one of
        the subclasses of ``Climatology`` if remapping will be performed.

    remappedDataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology.  Available only after
        calling the ``remap_and_write`` method.

    remappedFileName : str
        The name of the file where ``remappedDataSet`` (the climatology data
        set after remapping) should be stored.  This attribute should be set by
        one of the subclasses of ``Climatology`` if remapping will be
        performed.

    Examples
    --------
    Here is an example of how to use a ``Climatology`` object within a
    task. In the example, ``ds`` is an xarray data set with a ``Time``
    dimension.  The result of the example code is that the xarray data set
    ``climatology.dataSet`` is available for computation and plotting.

    >>> from ..shared.climatology import Climatology
    >>> climatology = Climatology(
                task=self,
                monthNames='ANN')
    >>> climatology.compute(ds=ds)
    >>> print climatology.dataSet

    Computing a monthly climatology would look like this:
    >>> from ..shared.climatology import Climatology
    >>> climatology = Climatology(task=self)
    >>> climatology.compute_monthly(ds=ds)
    >>> print climatology.dataSet

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, task, monthNames=None, comparisonGrid=None):
        """
        Create a new climatology.

        Parameters
        ----------
        task : ``AnalysisTask`` object
            An analysis task for which the climatology is needed. ``task`` is
            used to get config options, define the component name, map streams,
            namelists and variables, etc.

        monthNames : str, optional
            The months that make up the climatology, used in constructing the
            names of cache files.  If provided, ``monthNames`` should be one of
            the keys of ``monthDictionary`` in the ``constants`` module.

        comparisonGrid : {'latlon'}, optional
            The name of the comparison grid to use for remapping (if any).
            If ``comparisonGrid=None`` (the default), no remapping will be
            performed.

        Raises
        ------
        ValueError
            If comarisonGrid does not describe a known comparions grid

        Authors
        -------
        Xylar Asay-Davis
        """
        self.task = task
        self.config = task.config
        self.calendar = task.calendar
        self.monthNames = monthNames
        if monthNames is None:
            self.monthValues = None
        else:
            monthValues = constants.monthDictionary[monthNames]
            if isinstance(monthValues, numbers.Integral):
                self.monthValues = [monthValues]
            else:
                self.monthValues = sorted(monthValues)
        if comparisonGrid == 'latlon':
            self.comparisonDescriptor = \
                get_lat_lon_comparison_descriptor(task.config)
        elif comparisonGrid is None:
            self.comparisonDescriptor = None
        else:
            raise ValueError('Unknown comaprison grid type {}'.format(
                comparisonGrid))

    def create_remapper(self, mappingFileSection, mappingFileOption,
                        mappingFilePrefix, method):  # {{{
        """
        Creates an attribute ``remapper``, a ``Remapper`` object, that can be
        used to remap from source files or data sets to corresponding data sets
        on the comparison grid.

        This call requires that the ``sourceDescriptor`` attribute be assigned,
        which happens automatically for some of the subclasses of
        ``Climatology`` (e.g. ``MpasClimatology`` and
        ``ObservationClimatology``).

        If necessary, creates the mapping file containing weights and indices
        needed to perform remapping.

        Parameters
        ----------
        mappingFileSection, mappingFileOption : str
            Section and option in ``config`` where the name of the mapping file
            may be given, or where it will be stored if a new mapping file is
            created

        mappingFilePrefix : str
            A prefix to be prepended to the mapping file name

        method : {'bilinear', 'neareststod', 'conserve'}
            The method of interpolation used.

        Returns
        -------
        remapper : ``Remapper`` object
            A remapper between the source and comparison grids, also stored
            as a ``remapper`` attribute.

        Authors
        -------
        Xylar Asay-Davis
        """

        config = self.config

        if config.has_option(mappingFileSection, mappingFileOption):
            # a mapping file was supplied, so we'll use that name
            mappingFileName = config.get(mappingFileSection, mappingFileOption)
        else:
            if self._matches_comparison(self.sourceDescriptor,
                                        self.comparisonDescriptor):
                # no need to remap
                mappingFileName = None

            else:
                # we need to build the path to the mapping file and an
                # appropriate file name
                mappingSubdirectory = \
                    build_config_full_path(config, 'output',
                                           'mappingSubdirectory')

                make_directories(mappingSubdirectory)

                mappingFileName = '{}/{}_{}_to_{}_{}.nc'.format(
                    mappingSubdirectory, mappingFilePrefix,
                    self.sourceDescriptor.meshName,
                    self.comparisonDescriptor.meshName,
                    method)

                config.set(mappingFileSection, mappingFileOption,
                           mappingFileName)

        remapper = Remapper(self.sourceDescriptor, self.comparisonDescriptor,
                            mappingFileName)

        remapper.build_mapping_file(method=method)

        self.remapper = remapper
        return remapper  # }}}

    def compute(self, ds, monthValues=None, maskVaries=True):  # {{{
        """
        Compute a monthly, seasonal or annual climatology data set from a data
        set.  The mean is weighted but the number of days in each month of
        the data set, ignoring values masked out with NaNs.  If the month
        coordinate is not present, a data array ``month`` will be added based
        on ``Time`` and the provided calendar.

        Parameters
        ----------
        ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
            A data set with a ``Time`` coordinate expressed as days since
            0001-01-01 or ``month`` coordinate

        monthValues : int or array-like of ints, optional
            A single month or an array of months to be averaged together.
            If this option is not provided, the value of ``monthValues`` passed
            to ``__init__`` will be used.

        maskVaries: bool, optional
            If the mask (where variables in ``ds`` are ``NaN``) varies with
            time. If not, the weighted average does not need make extra effort
            to account for the mask.  Most MPAS fields will have masks that
            don't vary in time, whereas observations may sometimes be present
            only at some times and not at others, requiring
            ``maskVaries = True``.

        Returns
        -------
        dataSet : object of same type as ``ds``
            A data set without the ``'Time'`` coordinate containing the mean
            of ds over all months in monthValues, weighted by the number of
            days in each month.  Also stored as the ``dataSet`` attribute.

        Authors
        -------
        Xylar Asay-Davis
        """

        if monthValues is None:
            monthValues = self.monthValues

        ds = add_years_months_days_in_month(ds, self.calendar)

        mask = xr.zeros_like(ds.month, bool)

        for month in monthValues:
            mask = xr.ufuncs.logical_or(mask, ds.month == month)

        climatologyMonths = ds.where(mask, drop=True)

        self.dataSet = self._compute_masked_mean(climatologyMonths, maskVaries)

        return self.dataSet  # }}}

    def compute_monthly(self, ds, maskVaries=True):  # {{{
        """
        Compute monthly climatologies from a data set.  The mean is
        weighted by the number of days in each month of the data set,
        ignoring values masked out with NaNs.  If the month coordinate is
        not present, a data array ``month`` will be added based on
        ``Time``.

        Parameters
        ----------
        ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
            A data set with a ``Time`` coordinate expressed as days since
            0001-01-01 or ``month`` coordinate

        maskVaries: bool, optional
            If the mask (where variables in ``ds`` are ``NaN``) varies with
            time. If not, the weighted average does not need make extra
            effort to account for the mask.  Most MPAS fields will have
            masks that don't vary in time, whereas observations may
            sometimes be present only at some times and not at others,
            requiring ``maskVaries = True``.

        Returns
        -------
        dataSet : object of same type as ``ds``
            A data set with the same fields at ``ds`` but in which all data
            from each month has been averaged together to create a
            climatology for each month.  Also stored as the ``dataSet``
            attribute

        Authors
        -------
        Xylar Asay-Davis
        """

        def compute_one_month_climatology(ds):
            monthValues = list(ds.month.values)
            return self.compute(ds, monthValues, maskVaries)

        ds = add_years_months_days_in_month(ds, self.calendar)

        self.dataSet = \
            ds.groupby('month').apply(compute_one_month_climatology)

        return self.dataSet  # }}}

    def remap_and_write(self, useNcremap=None):  # {{{
        """
        Remap the climatology data set produced by a call to ``compute`` or
        ``cache``, write the result to an output file, and return the remapped
        data set.

        This call requires that ``climatologyFileName`` and
        ``remappedFileName`` attributes be assigned, which happens
        automatically for some of the subclasses of ``Climatology``
        (e.g. ``MpasClimatology`` and ``ObservationClimatology``).

        Note that the files named by attributes ``climatologyFileName`` and
        ``remappedFileName`` will be overwritten if they exist, so if
        this behavior is not desired, the calling code should skip this call
        if the files exist and simply load the contents of
        ``remappedFileName``.

        Parameters
        ----------
        useNcremap : bool, optional
            If present, overrides ``useNcremap`` from the config file.  This
            is useful for data sets that cannot be handle by ncremap (or are
            masked better with the "online" remapper).

        Returns
        -------
        remappedDataSet : ``xarray.DataSet`` or ``xarray.DataArray`` object
            A data set containing the remapped climatology.  Also stored in
            the ``remappedDataSet`` attribute.

        Authors
        -------
        Xylar Asay-Davis
        """
        if useNcremap is None:
            useNcremap = self.config.getboolean('climatology', 'useNcremap')

        if self.remapper.mappingFileName is None:
            # no remapping is needed
            self.remappedDataSet = self.dataSet
        else:
            if useNcremap:
                if not os.path.exists(self.climatologyFileName):
                    self.dataSet.to_netcdf(self.climatologyFileName)
                self.remapper.remap_file(inFileName=self.climatologyFileName,
                                         outFileName=self.remappedFileName,
                                         overwrite=True)
                self.remappedDataSet = xr.open_dataset(self.remappedFileName)
            else:
                renormalizationThreshold = self.config.getfloat(
                    'climatology', 'renormalizationThreshold')

                self.remappedDataSet = self.remapper.remap(
                    self.dataSet, renormalizationThreshold)
                self.remappedDataSet.to_netcdf(self.remappedFileName)
        return self.remappedDataSet  # }}}

    def _compute_masked_mean(self, ds, maskVaries):  # {{{
        '''
        Compute the time average of data set, masked out where the variables
        in ds are NaN and, if ``maskVaries == True``, weighting by the number
        of days used to compute each monthly mean time in ds.

        Authors
        -------
        Xylar Asay-Davis
        '''
        def ds_to_weights(ds):
            # make an identical data set to ds but replacing all data arrays
            # with notnull applied to that data array
            weights = ds.copy(deep=True)
            if isinstance(ds, xr.core.dataarray.DataArray):
                weights = ds.notnull()
            elif isinstance(ds, xr.core.dataset.Dataset):
                for var in ds.data_vars:
                    weights[var] = ds[var].notnull()
            else:
                raise TypeError('ds must be an instance of either '
                                'xarray.Dataset or xarray.DataArray.')

            return weights

        if maskVaries:
            dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time',
                                                      keep_attrs=True)

            weights = ds_to_weights(ds)

            weightSum = (weights * ds.daysInMonth).sum(dim='Time')

            timeMean = dsWeightedSum / weightSum.where(weightSum > 0.)
        else:
            days = ds.daysInMonth.sum(dim='Time')

            dsWeightedSum = (ds * ds.daysInMonth).sum(dim='Time',
                                                      keep_attrs=True)

            timeMean = dsWeightedSum / days.where(days > 0.)

        return timeMean  # }}}

    def _matches_comparison(self, obsDescriptor, comparisonDescriptor):  # {{{
        '''
        Determine if the two meshes are the same

        Authors
        -------
        Xylar Asay-Davis
        '''

        if isinstance(obsDescriptor, ProjectionGridDescriptor) and \
                isinstance(comparisonDescriptor, ProjectionGridDescriptor):
            # pretty hard to determine if projections are the same, so we'll
            # rely on the grid names
            match = \
                obsDescriptor.meshName == comparisonDescriptor.meshName and \
                len(obsDescriptor.x) == len(comparisonDescriptor.x) and \
                len(obsDescriptor.y) == len(comparisonDescriptor.y) and \
                numpy.all(numpy.isclose(obsDescriptor.x,
                                        comparisonDescriptor.x)) and \
                numpy.all(numpy.isclose(obsDescriptor.y,
                                        comparisonDescriptor.y))
        elif isinstance(obsDescriptor, LatLonGridDescriptor) and \
                isinstance(comparisonDescriptor, LatLonGridDescriptor):
            match = \
                ((('degree' in obsDescriptor.units and
                   'degree' in comparisonDescriptor.units) or
                  ('radian' in obsDescriptor.units and
                   'radian' in comparisonDescriptor.units)) and
                 len(obsDescriptor.lat) == len(comparisonDescriptor.lat) and
                 len(obsDescriptor.lon) == len(comparisonDescriptor.lon) and
                 numpy.all(numpy.isclose(obsDescriptor.lat,
                                         comparisonDescriptor.lat)) and
                 numpy.all(numpy.isclose(obsDescriptor.lon,
                                         comparisonDescriptor.lon)))
        else:
            match = False

        return match  # }}}

# }}}


class MpasClimatology(Climatology):  # {{{
    """
    A class for computing climatologies from an MPAS monthly mean data set.

    Attributes
    ----------
    fieldName : str
        The name of the field for which the climatology is being computed,
        used in constructing the names of cache files.

    streamName : str
        The name of the stream from which the climatology data set will
        be read, used to cache times and update their bounds in the
        configuration parser.

    startYear, endYear : int
        The start and end years of the climatology

    climatologyDirectory : str
        The directory where climatologies will be cached

    climatologyPrefix : str
        The prefix (including full path) on climatology cache files.  A call
        to the ``cache`` method may produce several intermediate cache files
        as well as a "final" cache file with the seasonal average between
        ``startYear`` and ``endYear``.

    climatologyFileName : str
        The full path to final climatology cache file where ``dataSet``
        (the climatology data set before remapping) should be stored.

    remappedFileName : str
        The name of the cache file where ``remappedDataSet`` (the climatology
        data set after remapping) should be stored.  This attribute will be
        set only if a comparison grid name is supplied.

    sourceDescriptor : ``MpasMeshDescriptor`` object
        A descriptor of the source MPAS mesh used for remapping.

    remapper : ``Remapper`` object
        A remapper between the source and comparison grids, created
        automatically if a comparison grid is supplied.

    dataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the climatology (before remapping, if
        applicable).  Available only after calling the ``compute``,
        ``compute_monthly``, or ``cache`` methods.

    remappedDataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology. Loaded automatically
        if a comparison grid is supplied and the cache file in
        ``remappedFileName`` already exists.  Otherwise, available only after
        calling the ``remap_and_write`` method.

    Examples
    --------
    Here is an example of how to use an ``MpasClimatology`` object within a
    task. In the example, ``self._open_mpas_dataset_part`` is a function with
    arguments ``inputFileNames``, ``startDate``, and ``endDate`` used to open
    and return part of the data set during caching.  The result of the example
    code is that the xarray data set ``climatology.remappedDataSet`` is
    available for computation and plotting, having either been read from a
    cache file or computed for the MPAS monthly time series data set.

    >>> from ..shared.climatology import MpasClimatology
    >>> restartFileName = self.runStreams.readpath('restart')[0]
    >>> climatology = MpasClimatology(
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
    >>> if climatology.remappedDataSet is None:
            climatology.cache(
                    openDataSetFunc=self._open_mpas_dataset_part,
                    printProgress=True)
            climatology.remap_and_write()
    >>> print climatology.remappedDataSet

    Here is another example, this time without remapping:

    >>> from ..shared.climatology import MpasClimatology
    >>> climatology = MpasClimatology(
                task=self,
                fieldName='sst',
                monthNames='ANN',
                streamName='timeSeriesStats')
    >>> climatology.cache(
                openDataSetFunc=self._open_mpas_dataset_part,
                printProgress=True)
    >>> print climatology.dataSet


    Authors
    -------
    Xylar Asay-Davis
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

        Parameters
        ----------
        task : ``AnalysisTask`` object
            An analysis task for which the climatology is needed. ``task`` is
            used to get config options, define the component name, map streams,
            namelists and variables, etc.

        fieldName : str
            The name of the field for which the climatology is being computed,
            used in constructing the names of cache files.

        monthNames : str
            The months that make up the climatology, used in constructing the
            names of cache files.

        streamName : str
            The name of the stream from which the climatology data set will
            be read, used to cache times and update their bounds in the
            configuration parser.

        The remaining parameters are all required if remapping is to be
        performed.

        meshFileName : str, optional
            The name of the MPAS mesh file, used to create a descriptor of the
            MPAS mesh for remapping.

        comparisonGrid : {'latlon'}, optional
            The name of the comparison grid to use for remapping (if any).

        mappingFileSection, mappingFileOption : str, optional
            Section and option in ``config`` where the name of the mapping file
            may be given, or where it will be stored if a new mapping file is
            created

        mappingFilePrefix : str, optional
            A prefix to be prepended to the mapping file name

        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used.

        Raises
        ------
        ValueError
            If comarisonGrid does not describe a known comparions grid

        Authors
        -------
        Xylar Asay-Davis
        """

        super(MpasClimatology, self).__init__(task, monthNames, comparisonGrid)

        self.fieldName = fieldName
        self.streamName = streamName
        climSection = 'climatology'
        task.update_start_end_date(section=climSection,
                                   streamName=streamName)

        self.startYear = self.config.getint(climSection, 'startYear')
        self.endYear = self.config.getint(climSection, 'endYear')

        mpasMeshName = self.config.get('input', 'mpasMeshName')

        self.climatologyDirectory = build_config_full_path(
            self.config, 'output', 'mpasClimatologySubdirectory')

        make_directories(self.climatologyDirectory)

        self.climatologyPrefix = \
            '{}/{}_{}_{}'.format(self.climatologyDirectory, fieldName,
                                 mpasMeshName, monthNames)

        yearString, fileSuffix = self._get_year_string(self.startYear,
                                                       self.endYear)
        self.climatologyFileName = \
            '{}_{}.nc'.format(self.climatologyPrefix, fileSuffix)

        if comparisonGrid is not None:
            remappedDirectory = build_config_full_path(
                self.config, 'output', 'mpasRemappedClimSubdirectory')

            make_directories(remappedDirectory)

            self.remappedFileName = '{}/{}_{}_to_{}_{}_{}.nc'.format(
                    remappedDirectory, fieldName, mpasMeshName,
                    self.comparisonDescriptor.meshName, monthNames, fileSuffix)

            if os.path.exists(self.remappedFileName):
                # no need to create the remapper, since the cached data set
                # alredy exists
                self.remappedDataSet = xr.open_dataset(self.remappedFileName)
            else:
                self.sourceDescriptor = \
                    MpasMeshDescriptor(fileName=meshFileName,
                                       meshName=mpasMeshName)

                self.create_remapper(mappingFileSection, mappingFileOption,
                                     mappingFilePrefix, method)

                self.remappedDataSet = None
        # }}}

    def cache(self, openDataSetFunc, printProgress=False):  # {{{
        '''
        Cache NetCDF files for each year of an annual climatology, and then use
        the cached files to compute a climatology for the full range of years.
        The start and end years of the climatology are taken from ``config``,
        and are updated in ``config`` if the data set ``ds`` doesn't contain
        this full range.

        Note: only works with climatologies where the mask (locations of
        ``NaN`` values) doesn't vary with time.

        Parameters
        ----------
        openDataSetFunc : function
            A function with arguments ``inputFileNames``, ``startDate``,
            and ``endDate`` used to open a portion of the data set from which
            the climatology is computed and cached.  Typically, the function
            makes a call to ``generalized_reader.open_multifile_dataset``,
            perhaps performing further manipulation of the resulting data set.

        printProgress: bool, optional
            Whether progress messages should be printed as the climatology is
            computed

        Returns
        -------
        dataSet : object of same type as ``ds``
            A data set without the ``'Time'`` coordinate containing the mean
            of ds over all months in monthValues, weighted by the number of
            days in each month.  Also avialable through the ``dataSet``
            attribute.

        Authors
        -------
        Xylar Asay-Davis
        '''
        cacheInfo = self._setup_climatology_caching(printProgress)

        # compute and store each cache file with interval yearsPerCacheFile
        self._cache_individual_climatologies(openDataSetFunc, cacheInfo,
                                             printProgress)

        # compute the aggregate climatology
        self.dataSet = self._cache_aggregated_climatology(cacheInfo,
                                                          printProgress)

        return self.dataSet  # }}}

    def _setup_climatology_caching(self, printProgress):  # {{{
        '''
        Determine which cache files already exist, which are incomplete and
        which years are present in each cache file (whether existing or to be
        created).

        Authors
        -------
        Xylar Asay-Davis
        '''
        yearsPerCacheFile = self.config.getint('climatology',
                                               'yearsPerCacheFile')
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

        if printProgress:
            print '   Computing and caching climatologies covering {}-year ' \
                  'spans...'.format(yearsPerCacheFile)

        cacheInfo = []

        firstFile = None
        lastFile = None
        # figure out which files to load and which years go in each file
        for firstYear in range(self.startYear, self.endYear+1,
                               yearsPerCacheFile):
            years = numpy.arange(firstYear, firstYear+yearsPerCacheFile)

            monthsIfDone = len(self.monthValues)*len(years)

            yearString, fileSuffix = self._get_year_string(years[0], years[-1])
            outputFileName = '{}_{}.nc'.format(self.climatologyPrefix,
                                               fileSuffix)

            done = False
            if os.path.exists(outputFileName):
                # already cached
                dsCached = None
                try:
                    dsCached = xr.open_dataset(outputFileName)
                except IOError:
                    # assuming the cache file is corrupt, so deleting it.
                    message = 'Deleting cache file {}, which appears to have' \
                              ' been corrupted.'.format(outputFileName)
                    warnings.warn(message)
                    os.remove(outputFileName)

                if ((dsCached is not None) and
                        (dsCached.attrs['totalMonths'] == monthsIfDone)):
                    # also complete, so we can move on
                    done = True
                if dsCached is not None:
                    dsCached.close()

            inputFileNames = []
            for year in years:
                for month in self.monthValues:
                    for fileName in timeCache:
                        timeCacheYears = timeCache[fileName]['years']
                        timeCacheMonths = timeCache[fileName]['months']
                        mask = numpy.logical_and(timeCacheYears == year,
                                                 timeCacheMonths == month)
                        if numpy.count_nonzero(mask) > 0:
                            inputFileNames.append(fileName)

            if len(inputFileNames) > 0:
                cacheDict = {'outputFileName': outputFileName,
                             'done': done,
                             'years': years,
                             'yearString': yearString,
                             'inputFileNames': inputFileNames}
                cacheInfo.append(cacheDict)
                if not done:
                    lastFile = inputFileNames[-1]
                    if firstFile is None:
                        firstFile = inputFileNames[0]

        if printProgress and firstFile is not None:
            print '\n  Caching data from files:\n' \
                  '    {} through\n    {}'.format(
                      os.path.basename(firstFile),
                      os.path.basename(lastFile))

        return cacheInfo  # }}}

    def _cache_individual_climatologies(self, openDataSetFunc, cacheInfo,
                                        printProgress):  # {{{
        '''
        Cache individual climatologies for later aggregation.

        Authors
        -------
        Xylar Asay-Davis
        '''

        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')

        for info in cacheInfo:
            if info['done']:
                continue
            outputFileName = info['outputFileName']
            yearString = info['yearString']
            years = info['years']
            inputFileNames = info['inputFileNames']
            ds = openDataSetFunc(inputFileNames, startDate, endDate)
            ds = add_years_months_days_in_month(ds, self.calendar)
            mask = numpy.zeros(ds.dims['Time'], bool)
            for year in years:
                mask = numpy.logical_or(mask, ds.year.values == year)
            ds.coords['cacheMask'] = ('Time', mask)
            dsYear = ds.where(ds.cacheMask, drop=True)

            if printProgress:
                print '     {}'.format(yearString)

            totalDays = dsYear.daysInMonth.sum(dim='Time').values

            monthCount = dsYear.dims['Time']

            climatology = self.compute(dsYear, maskVaries=False)

            climatology.attrs['totalDays'] = totalDays
            climatology.attrs['totalMonths'] = monthCount
            climatology.attrs['fingerprintClimo'] = fingerprint_generator()

            climatology.to_netcdf(outputFileName)
            climatology.close()

        # }}}

    def _cache_aggregated_climatology(self, cacheInfo, printProgress):  # {{{
        '''
        Cache aggregated climatology from individual climatologies.

        Authors
        -------
        Xylar Asay-Davis

        '''
        yearString, fileSuffix = self._get_year_string(self.startYear,
                                                       self.endYear)
        outputFileName = '{}_{}.nc'.format(self.climatologyPrefix, fileSuffix)

        done = False
        if len(cacheInfo) == 0:
            climatology = None
            done = True

        if os.path.exists(outputFileName):
            # already cached
            climatology = None
            try:
                climatology = xr.open_dataset(outputFileName)

            except IOError:
                # assuming the cache file is corrupt, so deleting it.
                message = 'Deleting cache file {}, which appears to have ' \
                          'been corrupted.'.format(outputFileName)
                warnings.warn(message)
                os.remove(outputFileName)

            if len(cacheInfo) == 1 and \
                    outputFileName == cacheInfo[0]['outputFileName']:
                # theres only one cache file and it already has the same name
                # as the aggregated file so no need to aggregate
                done = True

            elif climatology is not None:
                monthsIfDone = (self.endYear-self.startYear+1) * \
                    len(self.monthValues)
                if climatology.attrs['totalMonths'] == monthsIfDone:
                    # also complete, so we can move on
                    done = True
                else:
                    climatology.close()

        if not done:
            if printProgress:
                print '   Computing aggregated climatology ' \
                      '{}...'.format(yearString)

            first = True
            for info in cacheInfo:
                inputFileName = info['outputFileName']
                ds = xr.open_dataset(inputFileName)
                days = ds.attrs['totalDays']
                months = ds.attrs['totalMonths']
                if first:
                    totalDays = days
                    totalMonths = months
                    climatology = ds * days
                    first = False
                else:
                    totalDays += days
                    totalMonths += months
                    climatology = climatology + ds * days

                ds.close()
            climatology = climatology / totalDays

            climatology.attrs['totalDays'] = totalDays
            climatology.attrs['totalMonths'] = totalMonths
            climatology.attrs['fingerprintClimo'] = fingerprint_generator()

            climatology.to_netcdf(outputFileName)

        return climatology  # }}}

    def _get_year_string(self, startYear, endYear):  # {{{
        if startYear == endYear:
            yearString = '{:04d}'.format(startYear)
            fileSuffix = 'year{}'.format(yearString)
        else:
            yearString = '{:04d}-{:04d}'.format(startYear, endYear)
            fileSuffix = 'years{}'.format(yearString)

        return yearString, fileSuffix  # }}}

    # }}}


class ObservationClimatology(Climatology):  # {{{
    """
    A class for computing climatologies from an obsevational data set.

    Attributes
    ----------
    fieldName : str
        The name of the field for which the climatology is being computed,
        used in constructing the names of cache files.

    climatologyDirectory : str
        The directory where climatologies will be cached

    climatologyFileName : str
        The full path to final climatology cache file where ``dataSet``
        (the climatology data set before remapping) should be stored.

    remappedFileName : str
        The name of the cache file where ``remappedDataSet`` (the climatology
        data set after remapping) should be stored.  This attribute will be
        set only if a comparison grid name is supplied.

    sourceDescriptor : ``MeshDescriptor`` object
        The descriptor of the input data grid used for remapping.

    remapper : ``Remapper`` object
        A remapper between the source and comparison grids, created
        automatically if a comparison grid is supplied.

    dataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the climatology (before remapping, if
        applicable).  Available only after calling the ```compute`` or
        ``compute_monthly`` methods.

    remappedDataSet : ``xarray.Dataset`` or ``xarray.DataArray`` object
        A data set containing the remapped climatology. Loaded automatically
        if a comparison grid is supplied and the cache file in
        ``remappedFileName`` already exists.  Otherwise, available only after
        calling the ``remap_and_write`` method.

    Examples
    --------
    Here is an example of how to use an ``ObservationClimatology`` object
    within a task. In the example, ``obsDescriptor`` is a ``MeshDescriptor``
    object describing the input observations grid and ``dsObs`` is an xarray
    data set containing the monthly mean time series of the observations.  The
    result of the example code is that the xarray data set
    ``climatology.remappedDataSet`` is available for computation and plotting,
    having either been read from a cache file or computed for the monthly
    time series data set.

    >>> from ..shared.climatology import ObservationClimatology
    >>> climatology = ObservationClimatology(
                task=self,
                fieldName='sst',
                monthNames='ANN',
                obsGridDescriptor=obsDescriptor,
                comparisonGrid='latlon',
                mappingFileSection='oceanObservations',
                mappingFileOption='sstClimatologyMappingFile',
                mappingFilePrefix='map_obs_sst',
                method='bilinear')
    >>> if climatology.remappedDataSet is None:
            climatology.compute(ds=dsObs)
            climatology.remap_and_write()
    >>> print climatology.remappedDataSet


    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, task, fieldName, monthNames,
                 obsGridDescriptor, comparisonGrid=None,
                 mappingFileSection=None, mappingFileOption=None,
                 mappingFilePrefix=None, method=None):  # {{{
        """
        Create a new object for creating climatologies from observational data
        sets.

        Parameters
        ----------
        task : ``AnalysisTask`` object
            An analysis task for which the climatology is needed. ``task`` is
            used to get config options, define the component name, map streams,
            namelists and variables, etc.

        fieldName : str
            The name of the field for which the climatology is being computed,
            used in constructing the names of cache files.

        monthNames : str
            The months that make up the climatology, used in constructing the
            names of cache files.

        obsGridDescriptor : ``MeshDescriptor`` object
            The descriptor of the input data grid

        The remaining parameters are all required if remapping is to be
        performed.

        comparisonGrid : {'latlon'}, optional
            The name of the comparison grid to use for remapping (if any).

        mappingFileSection, mappingFileOption : str, optional
            Section and option in ``config`` where the name of the mapping file
            may be given, or where it will be stored if a new mapping file is
            created

        mappingFilePrefix : str, optional
            A prefix to be prepended to the mapping file name

        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used.

        Raises
        ------
        ValueError
            If comarisonGrid does not describe a known comparions grid

        Authors
        -------
        Xylar Asay-Davis
        """

        super(ObservationClimatology, self).__init__(task, monthNames,
                                                     comparisonGrid)

        self.fieldName = fieldName

        obsSection = '{}Observations'.format(task.componentName)

        climatologyDirectory = build_config_full_path(
            config=self.config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection=obsSection)

        self.sourceDescriptor = obsGridDescriptor

        obsGridName = self.sourceDescriptor.meshName
        comparisonGridName = self.comparisonDescriptor.meshName

        self.climatologyFileName = '{}/{}_{}_{}.nc'.format(
            climatologyDirectory, fieldName, obsGridName, monthNames)

        make_directories(climatologyDirectory)

        if comparisonGrid is not None:
            remappedDirectory = build_config_full_path(
                config=self.config, section='output',
                relativePathOption='remappedClimSubdirectory',
                relativePathSection=obsSection)

            if self._matches_comparison(self.sourceDescriptor,
                                        self.comparisonDescriptor):
                self.remappedFileName = self.climatologyFileName
            else:
                make_directories(remappedDirectory)

                self.remappedFileName = '{}/{}_{}_to_{}_{}.nc'.format(
                    remappedDirectory, fieldName, obsGridName,
                    comparisonGridName, monthNames)

            if os.path.exists(self.remappedFileName):
                # no need to create the remapper, since the cached data set
                # alredy exists
                self.remappedDataSet = xr.open_dataset(self.remappedFileName)
            else:

                self.create_remapper(mappingFileSection, mappingFileOption,
                                     mappingFilePrefix, method)

                self.remappedDataSet = None
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
