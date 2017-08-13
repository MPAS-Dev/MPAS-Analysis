import os
import os.path

import numpy.ma as ma
import numpy as np

import xarray as xr

from ..shared.constants import constants

from ..shared.climatology import get_lat_lon_comparison_descriptor, \
    get_remapper, get_mpas_climatology_dir_name, \
    get_observation_climatology_file_names, \
    compute_climatologies_with_ncclimo, \
    update_climatology_bounds_from_file_names, \
    remap_and_write_climatology, get_ncclimo_season_file_name
from ..shared.grid import MpasMeshDescriptor, LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path

from .sea_ice_analysis_task import SeaIceAnalysisTask


class ClimatologyMapSeaIce(SeaIceAnalysisTask):
    """
    General comparison of 2-d model fields against data.  Currently only
    supports sea ice concentration and sea ice thickness

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
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
        # first, call the constructor from the base class (SeaIceAnalysisTask)
        super(ClimatologyMapSeaIce, self).__init__(
            config=config,
            taskName='climatologyMapSeaIceConcThick',
            componentName='seaIce',
            tags=['climatology', 'horizontalMap'])

        # }}}

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
        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')
        self.inputFiles = \
            self.historyStreams.readpath(streamName,
                                         startDate=startDate,
                                         endDate=endDate,
                                         calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, startDate,
                                       endDate))

        update_climatology_bounds_from_file_names(self.inputFiles,
                                                  self.config)

        self.startYear = self.config.getint('climatology', 'startYear')
        self.endYear = self.config.getint('climatology', 'endYear')

        # }}}

    def run(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of sea-ice concentration and thickness " \
            "climatologies..."

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

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

        self._compute_seasonal_climatologies(
                seasons=['JFM', 'JAS', 'DJF', 'JJA', 'FM', 'ON'])

        self._compute_and_plot_concentration()

        self._compute_and_plot_thickness()   # }}}

    def _compute_seasonal_climatologies(self, seasons):  # {{{

        config = self.config

        mpasMeshName = self.mpasRemapper.sourceDescriptor.meshName
        comparisonGridName = self.mpasRemapper.destinationDescriptor.meshName

        overwriteMpasClimatology = config.getWithDefault(
            'climatology', 'overwriteMpasClimatology', False)

        self.climatologyDirectory = \
            get_mpas_climatology_dir_name(
                config=config,
                fieldName='seaIceAreaVol',
                mpasMeshName=mpasMeshName)

        (self.maskedClimatologyDirectory, self.remappedDirectory) = \
            get_mpas_climatology_dir_name(
                config=config,
                fieldName='seaIceAreaVol_masked',
                mpasMeshName=mpasMeshName,
                comparisonGridName=comparisonGridName)

        # the last climatology produced by NCO is always the annual, so if that
        # exists, others are also finished.
        climatologyFileName = \
            '{}/mpascice_ANN_climo.nc'.format(self.climatologyDirectory)
        if (overwriteMpasClimatology or
                not os.path.exists(climatologyFileName)):

            compute_climatologies_with_ncclimo(
                    config=config,
                    inDirectory=self.historyDirectory,
                    outDirectory=self.climatologyDirectory,
                    startYear=self.startYear,
                    endYear=self.endYear,
                    variableList=['timeMonthly_avg_iceAreaCell',
                                  'timeMonthly_avg_iceVolumeCell'],
                    modelName='mpascice',
                    seasons=seasons,
                    decemberMode='sdd')

        self._remap_seasonal_climatology(seasons)

        # }}}

    def _compute_and_plot_concentration(self):
        '''
        computes seasonal climatologies and plots model results, observations
        and biases in sea-ice concentration.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        '''

        print '  Make ice concentration plots...'

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')

        subtitle = 'Ice concentration'

        hemisphereSeasons = {'JFM': ('NH', 'Winter'),
                             'JAS': ('NH', 'Summer'),
                             'DJF': ('SH', 'Winter'),
                             'JJA': ('SH', 'Summer')}

        obsFileNames = {}
        remappedObsFileNames = {}
        obsRemappers = {}

        comparisonDescriptor = self.mpasRemapper.destinationDescriptor

        buildObsClimatologies = False
        for months in hemisphereSeasons:
            hemisphere, season = hemisphereSeasons[months]
            climFieldName = 'iceConcentration'
            for obsName in ['NASATeam', 'Bootstrap']:
                key = (months, obsName)
                obsFileName = build_config_full_path(
                    config, 'seaIceObservations',
                    'concentration{}{}_{}'.format(obsName, hemisphere, months))
                obsFieldName = '{}_{}_{}'.format(climFieldName, hemisphere,
                                                 obsName)

                obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                          latVarName='t_lat',
                                                          lonVarName='t_lon')
                obsRemapper = get_remapper(
                        config=config, sourceDescriptor=obsDescriptor,
                        comparisonDescriptor=comparisonDescriptor,
                        mappingFilePrefix='map_obs_seaIce',
                        method=config.get('seaIceObservations',
                                          'interpolationMethod'))
                obsRemappers[key] = obsRemapper

                if not os.path.isfile(obsFileName):
                    raise OSError('Obs file {} not found.'.format(
                        obsFileName))

                (climatologyFileName, remappedFileName) = \
                    get_observation_climatology_file_names(
                        config=config, fieldName=obsFieldName,
                        monthNames=months, componentName=self.componentName,
                        remapper=obsRemapper)

                obsFileNames[key] = obsFileName
                remappedObsFileNames[key] = remappedFileName

                if not os.path.exists(remappedFileName):
                    buildObsClimatologies = True

        for months in hemisphereSeasons:
            hemisphere, season = hemisphereSeasons[months]
            field = 'timeMonthly_avg_iceAreaCell'
            climFieldName = 'iceConcentration'

            remappedFileName = \
                get_ncclimo_season_file_name(self.remappedDirectory,
                                             'mpascice', months,
                                             self.startYear, self.endYear)

            remappedClimatology = xr.open_dataset(remappedFileName)

            iceConcentration = remappedClimatology[field].values
            lon = remappedClimatology['lon'].values
            lat = remappedClimatology['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            if hemisphere == 'NH':
                plotProjection = 'npstere'
            else:
                plotProjection = 'spstere'

            (colormapResult, colorbarLevelsResult) = setup_colormap(
                config,
                'climatologyMapSeaIceConcThick',
                suffix='ConcResult{}'.format(season))
            (colormapDifference, colorbarLevelsDifference) = setup_colormap(
                config,
                'climatologyMapSeaIceConcThick',
                suffix='ConcDifference{}'.format(season))

            referenceLongitude = config.getfloat(
                'climatologyMapSeaIceConcThick',
                'referenceLongitude{}'.format(hemisphere))
            minimumLatitude = config.getfloat(
                'climatologyMapSeaIceConcThick',
                'minimumLatitude{}'.format(hemisphere))

            # ice concentrations from NASATeam (or Bootstrap) algorithm
            for obsName in ['NASATeam', 'Bootstrap']:
                obsFieldName = 'AICE'

                key = (months, obsName)
                remappedFileName = remappedObsFileNames[key]

                if buildObsClimatologies:
                    obsFileName = obsFileNames[key]

                    seasonalClimatology = xr.open_dataset(obsFileName)

                    remappedClimatology = remap_and_write_climatology(
                            config, seasonalClimatology,  climatologyFileName,
                            remappedFileName, obsRemappers[key])

                obsIceConcentration = remappedClimatology[obsFieldName].values

                difference = iceConcentration - obsIceConcentration

                title = '{} ({}, years {:04d}-{:04d})'.format(
                    subtitle, months, self.startYear, self.endYear)
                fileout = '{}/iceconc{}{}_{}_{}_years{:04d}-{:04d}.png'.format(
                    self.plotsDirectory, obsName, hemisphere, mainRunName,
                    months, self.startYear, self.endYear)
                plot_polar_comparison(
                    config,
                    lonTarg,
                    latTarg,
                    iceConcentration,
                    obsIceConcentration,
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
                    obsTitle='Observations (SSM/I {})'.format(obsName),
                    diffTitle='Model-Observations',
                    cbarlabel='fraction')

    def _compute_and_plot_thickness(self):
        '''
        Computes seasonal climatologies and plots model results, observations
        and biases in sea-ice thickness.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        '''

        print '  Make ice thickness plots...'

        config = self.config

        subtitle = 'Ice thickness'

        plotsDirectory = build_config_full_path(config, 'output',
                                                'plotsSubdirectory')
        mainRunName = config.get('runs', 'mainRunName')

        obsFileNames = {}
        remappedObsFileNames = {}
        obsRemappers = {}

        comparisonDescriptor = self.mpasRemapper.destinationDescriptor

        # build a list of remapped observations files
        buildObsClimatologies = False
        for months in ['FM', 'ON']:
            climFieldName = 'iceThickness'
            for hemisphere in ['NH', 'SH']:
                key = (months, hemisphere)
                obsFileName = build_config_full_path(
                    config, 'seaIceObservations',
                    'thickness{}_{}'.format(hemisphere, months))
                if not os.path.isfile(obsFileName):
                    raise OSError('Obs file {} not found.'.format(
                        obsFileName))

                obsFieldName = '{}_{}'.format(climFieldName, hemisphere)
                obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                          latVarName='t_lat',
                                                          lonVarName='t_lon')
                obsRemapper = get_remapper(
                        config=config, sourceDescriptor=obsDescriptor,
                        comparisonDescriptor=comparisonDescriptor,
                        mappingFilePrefix='map_obs_seaIce',
                        method=config.get('seaIceObservations',
                                          'interpolationMethod'))
                obsRemappers[key] = obsRemapper

                (climatologyFileName, remappedFileName) = \
                    get_observation_climatology_file_names(
                        config=config, fieldName=obsFieldName,
                        monthNames=months, componentName=self.componentName,
                        remapper=obsRemapper)

                obsFileNames[key] = obsFileName
                remappedObsFileNames[key] = remappedFileName

                if not os.path.exists(remappedFileName):
                    buildObsClimatologies = True

        for months in ['FM', 'ON']:
            field = 'timeMonthly_avg_iceVolumeCell'
            climFieldName = 'iceThickness'

            remappedFileName = \
                get_ncclimo_season_file_name(self.remappedDirectory,
                                             'mpascice', months,
                                             self.startYear, self.endYear)

            remappedClimatology = xr.open_dataset(remappedFileName)

            iceThickness = remappedClimatology[field].values
            iceThickness = ma.masked_values(iceThickness, 0)
            lon = remappedClimatology['lon'].values
            lat = remappedClimatology['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            for hemisphere in ['NH', 'SH']:
                obsFieldName = 'HI'

                (colormapResult, colorbarLevelsResult) = setup_colormap(
                    config,
                    'climatologyMapSeaIceConcThick',
                    suffix='ThickResult{}'.format(hemisphere))
                (colormapDifference, colorbarLevelsDifference) = \
                    setup_colormap(config, 'climatologyMapSeaIceConcThick',
                                   suffix='ThickDifference{}'.format(
                                       hemisphere))

                referenceLongitude = config.getfloat(
                    'climatologyMapSeaIceConcThick',
                    'referenceLongitude{}'.format(hemisphere))
                minimumLatitude = config.getfloat(
                    'climatologyMapSeaIceConcThick',
                    'minimumLatitude{}'.format(hemisphere))

                # now the observations
                key = (months, hemisphere)
                remappedFileName = remappedObsFileNames[key]

                if buildObsClimatologies:
                    obsFileName = obsFileNames[key]

                    seasonalClimatology = xr.open_dataset(obsFileName)

                    remappedClimatology = remap_and_write_climatology(
                            config, seasonalClimatology, climatologyFileName,
                            remappedFileName, obsRemappers[key])

                obsIceThickness = remappedClimatology[obsFieldName].values

                # Mask thickness fields
                obsIceThickness = ma.masked_values(obsIceThickness, 0)
                if hemisphere == 'NH':
                    # Obs thickness should be nan above 86 (ICESat data)
                    obsIceThickness[latTarg > 86] = ma.masked
                    plotProjection = 'npstere'
                else:
                    plotProjection = 'spstere'

                difference = iceThickness - obsIceThickness

                title = '{} ({}, years {:04d}-{:04d})'.format(subtitle, months,
                                                              self.startYear,
                                                              self.endYear)
                fileout = '{}/icethick{}_{}_{}_years{:04d}-{:04d}.png'.format(
                    plotsDirectory, hemisphere, mainRunName, months,
                    self.startYear, self.endYear)
                plot_polar_comparison(
                    config,
                    lonTarg,
                    latTarg,
                    iceThickness,
                    obsIceThickness,
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
                    obsTitle='Observations (ICESat)',
                    diffTitle='Model-Observations',
                    cbarlabel='m')

        # }}}

    def _remap_seasonal_climatology(self, seasons):  # {{{

        modelName = 'mpascice'

        for season in seasons:

            overwriteMpasClimatology = self.config.getWithDefault(
                'climatology', 'overwriteMpasClimatology', False)

            monthValues = constants.monthDictionary[season]

            # interpolate the model results
            climatologyFileName = \
                get_ncclimo_season_file_name(self.climatologyDirectory,
                                             modelName, season,
                                             self.startYear, self.endYear)

            maskedClimatologyFileName = \
                get_ncclimo_season_file_name(self.maskedClimatologyDirectory,
                                             modelName, season,
                                             self.startYear, self.endYear)

            remappedFileName = \
                get_ncclimo_season_file_name(self.remappedDirectory,
                                             modelName, season,
                                             self.startYear, self.endYear)

            if (overwriteMpasClimatology or
                    not os.path.exists(maskedClimatologyFileName)):
                # slice the data set and set _FillValue (happens automatically)
                print climatologyFileName
                climatology = xr.open_dataset(climatologyFileName)
                iselValues = {'Time': 0}
                # select only Time=0
                climatology = climatology.isel(**iselValues)

                climatology.to_netcdf(maskedClimatologyFileName)

            if (overwriteMpasClimatology or
                    not os.path.exists(remappedFileName)):
                self.mpasRemapper.remap_file(
                        inFileName=maskedClimatologyFileName,
                        outFileName=remappedFileName,
                        overwrite=True)
        # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
