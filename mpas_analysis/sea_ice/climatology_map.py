import os
import os.path

import numpy.ma as ma
import numpy as np

import xarray as xr

from ..shared.constants import constants

from ..shared.climatology import get_lat_lon_comparison_descriptor, \
    get_remapper, get_mpas_climatology_file_names, \
    get_observation_climatology_file_names, \
    cache_climatologies, \
    update_climatology_bounds_from_file_names, \
    remap_and_write_climatology
from ..shared.grid import MpasMeshDescriptor, LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path
from ..shared.io import write_netcdf


from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from .sea_ice_analysis_task import SeaIceAnalysisTask


class ClimatologyMapSeaIce(SeaIceAnalysisTask):
    """
    General comparison of 2-d model fields against data.  Currently only
    supports sea ice concentration and sea ice thickness

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
    """

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (SeaIceAnalysisTask),
        # which will perform some common setup
        super(ClimatologyMapSeaIce, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # get a list of timeSeriesStatsMonthly output files from the streams
        # file, reading only those that are between the start and end dates
        streamName = self.historyStreams.find_stream(
            self.streamMap['timeSeriesStats'])
        self.startDate = self.config.get('climatology', 'startDate')
        self.endDate = self.config.get('climatology', 'endDate')
        self.inputFiles = \
            self.historyStreams.readpath(streamName,
                                         startDate=self.startDate,
                                         endDate=self.endDate,
                                         calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDate,
                                       self.endDate))

        changed, self.startYear, self.endYear, self.startDate, self.endDate = \
            update_climatology_bounds_from_file_names(self.inputFiles,
                                                      self.config)
        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of {} climatologies...".format(
            self.fieldNameInTitle)

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))
        # Load data

        print '  Load sea-ice data...'
        varList = [self.mpasFieldName]

        self.ds = open_multifile_dataset(
                fileNames=self.inputFiles,
                calendar=self.calendar,
                config=self.config,
                simulationStartTime=self.simulationStartTime,
                timeVariableName='Time',
                variableList=varList,
                iselValues=self.iselValues,
                variableMap=self.variableMap,
                startDate=self.startDate,
                endDate=self.endDate)

        # Compute climatologies (first motnhly and then seasonally)
        print '  Compute seasonal climatologies...'

        mpasDescriptor = MpasMeshDescriptor(
            self.restartFileName,
            meshName=self.config.get('input', 'mpasMeshName'))

        comparisonDescriptor = get_lat_lon_comparison_descriptor(self.config)

        parallel = self.config.getint('execute', 'parallelTaskCount') > 1
        if parallel:
            # avoid writing the same mapping file from multiple processes
            mappingFilePrefix = 'map_{}'.format(self.taskName)
        else:
            mappingFilePrefix = 'map'

        self.mpasRemapper = get_remapper(
            config=self.config, sourceDescriptor=mpasDescriptor,
            comparisonDescriptor=comparisonDescriptor,
            mappingFilePrefix=mappingFilePrefix,
            method=self.config.get('climatology', 'mpasInterpolationMethod'))

        self._compute_and_plot()  # }}}

    def _compute_and_plot(self):  # {{{
        '''
        computes seasonal climatologies and plots model results, observations
        and biases.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        '''

        print '  Make ice concentration plots...'

        config = self.config
        calendar = self.calendar
        ds = self.ds

        mainRunName = config.get('runs', 'mainRunName')
        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        comparisonDescriptor = self.mpasRemapper.destinationDescriptor

        hemisphere = self.hemisphere
        sectionName = self.sectionName
        vertical = config.getboolean(sectionName, 'vertical')
        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        for info in self.obsAndPlotInfo:
            season = info['season']
            monthValues = constants.monthDictionary[season]

            (colormapResult, colorbarLevelsResult) = setup_colormap(
                config, sectionName, suffix='Result')
            (colormapDifference, colorbarLevelsDifference) = setup_colormap(
                config, sectionName, suffix='Difference')

            referenceLongitude = config.getfloat(sectionName,
                                                 'referenceLongitude')
            minimumLatitude = config.getfloat(sectionName,
                                              'minimumLatitude')

            fieldName = '{}{}'.format(self.mpasFieldName, hemisphere)
            # interpolate the model results
            mpasMeshName = self.mpasRemapper.sourceDescriptor.meshName
            comparisonGridName = \
                self.mpasRemapper.destinationDescriptor.meshName
            (climatologyFileName, climatologyPrefix, remappedFileName) = \
                get_mpas_climatology_file_names(
                    config=config,
                    fieldName=fieldName,
                    monthNames=season,
                    mpasMeshName=mpasMeshName,
                    comparisonGridName=comparisonGridName)

            if not os.path.exists(climatologyFileName):
                seasonalClimatology = cache_climatologies(
                    ds, monthValues, config, climatologyPrefix, calendar,
                    printProgress=True)
                if seasonalClimatology is None:
                    # apparently, there was no data available to create the
                    # climatology
                    print 'Warning: no data to create sea ice thickness ' \
                          'climatology for {}'.format(season)
                    continue

                remappedClimatology = remap_and_write_climatology(
                    config, seasonalClimatology, climatologyFileName,
                    remappedFileName, self.mpasRemapper)

            else:

                remappedClimatology = xr.open_dataset(remappedFileName)

            modelOutput = remappedClimatology[self.mpasFieldName].values
            if self.maskValue is not None:
                modelOutput = ma.masked_values(modelOutput, self.maskValue)
            lon = remappedClimatology['lon'].values
            lat = remappedClimatology['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            obsFileName = info['obsFileName']

            obsDescriptor = LatLonGridDescriptor()
            obsDescriptor.read(fileName=obsFileName,
                               latVarName='t_lat',
                               lonVarName='t_lon')

            if not os.path.isfile(obsFileName):
                raise OSError('Obs file {} not found.'.format(
                    obsFileName))

            fieldName = '{}{}'.format(self.obsFieldName, hemisphere)
            obsRemapper = get_remapper(
                    config=config,
                    sourceDescriptor=obsDescriptor,
                    comparisonDescriptor=comparisonDescriptor,
                    mappingFilePrefix='map_obs_{}'.format(fieldName),
                    method=config.get('seaIceObservations',
                                      'interpolationMethod'))

            (obsClimatologyFileName, obsRemappedFileName) = \
                get_observation_climatology_file_names(
                    config=config, fieldName=fieldName,
                    monthNames=season, componentName=self.componentName,
                    remapper=obsRemapper)

            if not os.path.exists(obsRemappedFileName):

                # load the observations the first time
                seasonalClimatology = self._build_observational_dataset(
                        obsFileName)
                write_netcdf(seasonalClimatology, obsClimatologyFileName)

                if obsRemapper is None:
                    # no need to remap because the observations are on the
                    # comparison grid already
                    remappedClimatology = seasonalClimatology
                else:
                    remappedClimatology = \
                        remap_and_write_climatology(
                            config, seasonalClimatology,
                            obsClimatologyFileName,
                            obsRemappedFileName, obsRemapper)

            else:

                remappedClimatology = xr.open_dataset(obsRemappedFileName)

            observations = remappedClimatology[self.obsFieldName].values
            if self.maskValue is not None:
                observations = ma.masked_values(observations, self.maskValue)

            difference = modelOutput - observations

            startYear = self.startYear
            endYear = self.endYear
            outFileLabel = info['outFileLabel']
            observationTitleLabel = info['observationTitleLabel']

            title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, startYear, endYear)
            fileout = '{}/{}_{}_{}_years{:04d}-{:04d}.png'.format(
                self.plotsDirectory, outFileLabel,
                mainRunName, season, startYear, endYear)
            plot_polar_comparison(
                config,
                lonTarg,
                latTarg,
                modelOutput,
                observations,
                difference,
                colormapResult,
                colorbarLevelsResult,
                colormapDifference,
                colorbarLevelsDifference,
                title=title,
                fileout=fileout,
                plotProjection=plotProjection,
                latmin=minimumLatitude,
                lon0=referenceLongitude,
                modelTitle=mainRunName,
                obsTitle=observationTitleLabel,
                diffTitle='Model-Observations',
                cbarlabel=self.unitsLabel,
                vertical=vertical)
        # }}}
    # }}}


class ClimatologyMapSeaIceConc(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice concentration against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, hemisphere):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        hemisphere : ['NH', 'SH']
            The hemisphere to plot

        Authors
        -------
        Xylar Asay-Davis
        """
        taskName = 'climatologyMapSeaIceConc{}'.format(hemisphere)

        self.sectionName = taskName

        self.fieldName = 'seaIceConc'
        self.hemisphere = hemisphere
        self.fieldNameInTitle = 'Sea ice concentration'
        self.seasons = config.getExpression(self.sectionName, 'seasons')
        self.observationPrefixes = config.getExpression(self.sectionName,
                                                        'observationPrefixes')

        self.mpasFieldName = 'iceAreaCell'
        self.iselValues = None

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(config=config,
                                                       taskName=taskName,
                                                       componentName='seaIce',
                                                       tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceConc, self).setup_and_check()

        hemisphere = self.hemisphere

        self.obsAndPlotInfo = []
        for prefix in self.observationPrefixes:
            for season in self.seasons:
                localDict = {}
                localDict['season'] = season
                localDict['obsFileName'] = \
                    build_config_full_path(
                        self.config, 'seaIceObservations',
                        'concentration{}{}_{}'.format(prefix,
                                                      hemisphere,
                                                      season))
                localDict['observationTitleLabel'] = \
                    'Observations (SSM/I {})'.format(prefix)
                localDict['outFileLabel'] = \
                    'iceconc{}{}'.format(prefix, hemisphere)

                self.obsAndPlotInfo.append(localDict)

        self.obsFieldName = 'AICE'

        self.unitsLabel = 'fraction'

        self.maskValue = None

        # }}}

    def _build_observational_dataset(self, obsFileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(obsFileName)

        return dsObs  # }}}

    # }}}


class ClimatologyMapSeaIceThick(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice thickness against observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, hemisphere):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        hemisphere : ['NH', 'SH']
            The hemisphere to plot

        Authors
        -------
        Xylar Asay-Davis
        """
        taskName = 'climatologyMapSeaIceThick{}'.format(hemisphere)

        self.sectionName = taskName

        self.fieldName = 'seaIceThick'
        self.hemisphere = hemisphere
        self.fieldNameInTitle = 'Sea ice thickness'
        self.seasons = config.getExpression(self.sectionName, 'seasons')

        self.mpasFieldName = 'iceVolumeCell'
        self.iselValues = None

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(config=config,
                                                        taskName=taskName,
                                                        componentName='seaIce',
                                                        tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceThick, self).setup_and_check()

        hemisphere = self.hemisphere

        self.obsAndPlotInfo = []
        for season in self.seasons:
            localDict = {}
            localDict['season'] = season
            localDict['obsFileName'] = \
                build_config_full_path(
                    self.config, 'seaIceObservations',
                    'thickness{}_{}'.format(hemisphere, season))
            localDict['observationTitleLabel'] = 'Observations (ICESat)'
            localDict['outFileLabel'] = 'icethick{}'.format(hemisphere)

            self.obsAndPlotInfo.append(localDict)

        self.obsFieldName = 'HI'

        self.unitsLabel = 'm'

        self.maskValue = 0

        # }}}

    def _build_observational_dataset(self, obsFileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(obsFileName)

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
