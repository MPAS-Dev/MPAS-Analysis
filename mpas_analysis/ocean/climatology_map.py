"""
An analysis tasks for comparison of 2D model fields against observations.
Currently supports sea surface temperature (sst), sea surface salinity
(sss) and mixed layer depth (mld)

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
"""

import xarray as xr
import datetime
import numpy as np

from ..shared.analysis_task import AnalysisTask

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap
from ..shared.constants import constants

from ..shared.io import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.grid import LatLonGridDescriptor
from ..shared.climatology import MpasClimatology, ObservationClimatology

from ..shared.mpas_xarray import mpas_xarray


class ClimatologyMapOcean(AnalysisTask):  # {{{
    """
    An analysis task for comparison of 2D model fields against observations.

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, config, taskName, componentName, tags=[]):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        taskName :  str
            The name of the task, typically the same as the class name except
            starting with lowercase (e.g. 'myTask' for class 'MyTask')

        componentName :  {'ocean', 'seaIce'}
            The name of the component (same as the folder where the task
            resides)

        tags :  list of str, optional
            Tags used to describe the task (e.g. 'timeSeries', 'climatology',
            horizontalMap', 'index', 'transect').  These are used to determine
            which tasks are generated (e.g. 'all_transect' or 'no_climatology'
            in the 'generate' flags)

        Authors
        -------
        Xylar Asay-Davis
        """
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapOcean, self).__init__(
                config, taskName, componentName, tags,
                prerequisiteTasks=['cacheOceanTimeSeriesStatsTimes'])

        # by default, we don't override "useNcremap" from the config file.
        # However, some child classes will set this explicitly, since results
        # are better for remapping observations with the masking provided in
        # MPAS-Analysis directly
        self.useNcremapObs = None

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapOcean, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        self.simulationStartTime = get_simulation_start_time(self.runStreams)

        # }}}

    def run(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

        Last Modified
        -------------
        03/16/2017
        """

        print ""
        print "Plotting 2-d maps of {} climatologies...".format(
            self.fieldNameInTitle)

        # get local versions of member variables for convenience
        config = self.config
        fieldName = self.fieldName

        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for ocn_modelvsobs calculation')

        outputTimes = config.getExpression(self.taskName, 'comparisonTimes')

        parallel = self.config.getint('execute', 'parallelTaskCount') > 1
        if parallel:
            # avoid writing the same mapping file from multiple processes
            mappingFilePrefix = 'map_{}'.format(self.taskName)
        else:
            mappingFilePrefix = 'map'

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, self.taskName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, self.taskName, suffix='Difference')

        # we don't have any way to know the observation lat/lon and grid
        # without this, so we have to read the observational data set even
        # if it's already been remapped.
        dsObs = self._build_observational_dataset()
        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(ds=dsObs)

        # Interpolate and compute biases
        for months in outputTimes:
            monthValues = constants.monthDictionary[months]

            mpasClimatology = MpasClimatology(
                task=self,
                fieldName=fieldName,
                monthNames=months,
                streamName='timeSeriesStats',
                meshFileName=restartFileName,
                comparisonGrid='latlon',
                mappingFileSection='climatology',
                mappingFileOption='mpasMappingFile',
                mappingFilePrefix=mappingFilePrefix,
                method=config.get('climatology', 'mpasInterpolationMethod'))

            if mpasClimatology.remappedDataSet is None:
                # the remapped climatology hasn't been cached yet
                mpasClimatology.cache(
                    openDataSetFunc=self._open_mpas_dataset_part,
                    printProgress=True)
                mpasClimatology.remap_and_write()

            modelOutput = \
                mpasClimatology.remappedDataSet[self.mpasFieldName].values
            lon = mpasClimatology.remappedDataSet['lon'].values
            lat = mpasClimatology.remappedDataSet['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            obsClimatology = \
                ObservationClimatology(
                    task=self,
                    fieldName=self.obsFieldName,
                    monthNames=months,
                    obsGridDescriptor=obsDescriptor,
                    comparisonGrid='latlon',
                    mappingFileSection='oceanObservations',
                    mappingFileOption='{}ClimatologyMappingFile'.format(
                        fieldName),
                    mappingFilePrefix='map_obs_{}'.format(fieldName),
                    method=config.get('oceanObservations',
                                      'interpolationMethod'))

            if obsClimatology.remappedDataSet is None:
                # the remapped climatology hasn't been cached yet
                obsClimatology.compute(ds=dsObs, monthValues=monthValues,
                                       maskVaries=True)
                obsClimatology.remap_and_write(useNcremap=self.useNcremapObs)

            observations = \
                obsClimatology.remappedDataSet[self.obsFieldName].values

            bias = modelOutput - observations

            startYear = mpasClimatology.startYear
            endYear = mpasClimatology.endYear

            mainRunName = config.get('runs', 'mainRunName')
            outFileName = '{}/{}_{}_{}_years{:04d}-{:04d}.png'.format(
                    self.plotsDirectory, self.outFileLabel, mainRunName,
                    months, startYear, endYear)
            title = '{} ({}, years {:04d}-{:04d})'.format(
                    self.fieldNameInTitle, months, startYear, endYear)
            plot_global_comparison(config,
                                   lonTarg,
                                   latTarg,
                                   modelOutput,
                                   observations,
                                   bias,
                                   colormapResult,
                                   colorbarLevelsResult,
                                   colormapDifference,
                                   colorbarLevelsDifference,
                                   fileout=outFileName,
                                   title=title,
                                   modelTitle='{}'.format(mainRunName),
                                   obsTitle=self.observationTitleLabel,
                                   diffTitle='Model-Observations',
                                   cbarlabel=self.unitsLabel)

        # }}}

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

        Authors
        -------
        Xylar Asay-Davis
        """
        varList = [self.mpasFieldName]

        ds = open_multifile_dataset(
            fileNames=inputFileNames,
            calendar=self.calendar,
            config=self.config,
            simulationStartTime=self.simulationStartTime,
            timeVariableName='Time',
            variableList=varList,
            iselValues=self.iselValues,
            variableMap=self.variableMap,
            startDate=startDate,
            endDate=endDate)
        return ds  # }}}

    # }}}


class ClimatologyMapSST(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface temperature (sst) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sst'
        self.fieldNameInTitle = 'SST'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSST, self).__init__(
            config=config,
            taskName='climatologyMapSST',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSST, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(
                observationsDirectory)

        self.iselValues = {'nVertLevels': 0}
        self.mpasFieldName = 'temperature'

        self.obsFieldName = 'SST'

        climStartYear = self.config.getint('oceanObservations',
                                           'sstClimatologyStartYear')
        climEndYear = self.config.getint('oceanObservations',
                                         'sstClimatologyEndYear')

        if climStartYear < 1925:
            period = 'pre-industrial'
        else:
            period = 'present-day'

        # Set appropriate figure labels for SST
        self.observationTitleLabel = \
            'Observations (Hadley/OI, {} {:04d}-{:04d})'.format(period,
                                                                climStartYear,
                                                                climEndYear)
        self.outFileLabel = 'sstHADOI'
        self.unitsLabel = r'$^o$C'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        climStartYear = self.config.getint('oceanObservations',
                                           'sstClimatologyStartYear')
        climEndYear = self.config.getint('oceanObservations',
                                         'sstClimatologyEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_mfdataset(self.obsFileName)
        dsObs.rename({'time': 'Time'}, inplace=True)
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs  # }}}

    # }}}


class ClimatologyMapSSS(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface salinity (sss) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sss'
        self.fieldNameInTitle = 'SSS'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSSS, self).__init__(
            config=config,
            taskName='climatologyMapSSS',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSSS, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/Aquarius_V3_SSS_Monthly.nc'.format(
                observationsDirectory)

        self.iselValues = {'nVertLevels': 0}
        self.mpasFieldName = 'salinity'

        self.obsFieldName = 'SSS'

        self.observationTitleLabel = 'Observations (Aquarius, 2011-2014)'
        self.outFileLabel = 'sssAquarius'
        self.unitsLabel = 'PSU'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        dsObs = xr.open_mfdataset(self.obsFileName)
        dsObs.rename({'time': 'Time'}, inplace=True)
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs  # }}}

    # }}}


class ClimatologyMapMLD(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of mixed layer depth (mld) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """

        # We want to use the "online" remapper because it does a better job of
        # masking missing values for MLD.
        self.useNcremapObs = False

        self.fieldName = 'mld'
        self.fieldNameInTitle = 'MLD'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapMLD, self).__init__(
            config=config,
            taskName='climatologyMapMLD',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapMLD, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/holtetalley_mld_climatology.nc'.format(
                observationsDirectory)

        self.iselValues = None
        self.mpasFieldName = 'mld'

        self.obsFieldName = 'mld_dt_mean'

        # Set appropriate MLD figure labels
        self.observationTitleLabel = \
            'Observations (HolteTalley density threshold MLD)'
        self.outFileLabel = 'mldHolteTalleyARGO'
        self.unitsLabel = 'm'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        # Load MLD observational data
        dsObs = xr.open_mfdataset(self.obsFileName)

        # Increment month value to be consistent with the model output
        dsObs.iMONTH.values += 1
        # Rename the dimensions to be consistent with other obs. data sets
        dsObs.rename({'month': 'calmonth', 'lat': 'latCoord',
                      'lon': 'lonCoord'}, inplace=True)
        dsObs.rename({'iMONTH': 'Time', 'iLAT': 'lat', 'iLON': 'lon'},
                     inplace=True)

        # set the coordinates now that the dimensions have the same names
        dsObs.coords['lat'] = dsObs['latCoord']
        dsObs.coords['lon'] = dsObs['lonCoord']
        dsObs.coords['Time'] = dsObs['calmonth']
        dsObs.coords['month'] = ('Time', np.array(dsObs['calmonth'], int))

        # no meaningful year since this is already a climatology
        dsObs.coords['year'] = ('Time', np.ones(dsObs.dims['Time'], int))

        dsObs = mpas_xarray.subset_variables(dsObs, [self.obsFieldName,
                                                     'month'])

        return dsObs  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
