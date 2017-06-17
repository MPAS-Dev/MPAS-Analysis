import os
import os.path

import numpy.ma as ma
import numpy as np

import warnings
import xarray as xr

from ..shared.constants import constants

from ..shared.climatology import MpasClimatology, ObservationClimatology
from ..shared.grid import LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from .sea_ice_analysis_task import SeaIceAnalysisTask


class ClimatologyMapSeaIce(SeaIceAnalysisTask):  # {{{
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
            tags=['climatology', 'horizontalMap'],
            prerequisiteTasks=['cacheSeaIceTimeSeriesStatsTimes'])

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

        return  # }}}

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

        parallel = self.config.getint('execute', 'parallelTaskCount') > 1
        if parallel:
            # avoid writing the same mapping file from multiple processes
            self.mappingFilePrefix = 'map_{}'.format(self.taskName)
        else:
            self.mappingFilePrefix = 'map'

        self._compute_and_plot_concentration()

        self._compute_and_plot_thickness()   # }}}

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

        for months in hemisphereSeasons:
            hemisphere, season = hemisphereSeasons[months]
            # used in _open_sea_ice_part
            self.fieldName = 'iceAreaCell'
            climFieldName = 'iceConcentration'

            mpasClimatology = MpasClimatology(
                task=self,
                fieldName=climFieldName,
                monthNames=months,
                streamName='timeSeriesStats',
                meshFileName=self.restartFileName,
                comparisonGrid='latlon',
                mappingFileSection='climatology',
                mappingFileOption='mpasMappingFile',
                mappingFilePrefix=self.mappingFilePrefix,
                method=config.get('climatology', 'mpasInterpolationMethod'))

            if mpasClimatology.remappedDataSet is None:
                mpasClimatology.cache(openDataSetFunc=self._open_sea_ice_part,
                                      printProgress=True)
                if mpasClimatology.dataSet is None:
                    # apparently, there was no data available to create the
                    # climatology
                    warnings.warn('no data to create sea ice concentration '
                                  'climatology for {}'.format(months))
                    continue

                mpasClimatology.remap_and_write()

            iceConcentration = \
                mpasClimatology.remappedDataSet[self.fieldName].values
            lon = mpasClimatology.remappedDataSet['lon'].values
            lat = mpasClimatology.remappedDataSet['lat'].values

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
                obsFileName = build_config_full_path(
                    config, 'seaIceObservations',
                    'concentration{}{}_{}'.format(obsName, hemisphere, months))
                if not os.path.isfile(obsFileName):
                    raise OSError('Obs file {} not found.'.format(
                        obsFileName))

                obsFieldName = '{}_{}_{}'.format(climFieldName, hemisphere,
                                                 obsName)
                obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                          latVarName='t_lat',
                                                          lonVarName='t_lon')

                obsClimatology = \
                    ObservationClimatology(
                        task=self,
                        fieldName=obsFieldName,
                        monthNames=months,
                        obsGridDescriptor=obsDescriptor,
                        comparisonGrid='latlon',
                        mappingFileSection='seaIceObservations',
                        mappingFileOption='seaIceClimatologyMappingFile',
                        mappingFilePrefix='map_obs_seaIce',
                        method=config.get('seaIceObservations',
                                          'interpolationMethod'))

                if obsClimatology.remappedDataSet is None:
                    obsClimatology.dataSet = xr.open_dataset(obsFileName)
                    obsClimatology.remap_and_write()

                obsIceConcentration = \
                    obsClimatology.remappedDataSet['AICE'].values

                difference = iceConcentration - obsIceConcentration

                title = '{} ({}, years {:04d}-{:04d})'.format(
                    subtitle, months, mpasClimatology.startYear,
                    mpasClimatology.endYear)
                fileout = '{}/iceconc{}{}_{}_{}_years{:04d}-{:04d}.png'.format(
                    self.plotsDirectory, obsName, hemisphere, mainRunName,
                    months, mpasClimatology.startYear, mpasClimatology.endYear)
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

        for months in ['FM', 'ON']:
            self.fieldName = 'iceVolumeCell'
            climFieldName = 'iceThickness'

            mpasClimatology = MpasClimatology(
                task=self,
                fieldName=climFieldName,
                monthNames=months,
                streamName='timeSeriesStats',
                meshFileName=self.restartFileName,
                comparisonGrid='latlon',
                mappingFileSection='climatology',
                mappingFileOption='mpasMappingFile',
                mappingFilePrefix=self.mappingFilePrefix,
                method=config.get('climatology', 'mpasInterpolationMethod'))

            if mpasClimatology.remappedDataSet is None:
                mpasClimatology.cache(openDataSetFunc=self._open_sea_ice_part,
                                      printProgress=True)
                if mpasClimatology.dataSet is None:
                    # apparently, there was no data available to create the
                    # climatology
                    warnings.warn('no data to create sea ice concentration '
                                  'climatology for {}'.format(months))
                    continue

                mpasClimatology.remap_and_write()

            iceThickness = \
                mpasClimatology.remappedDataSet[self.fieldName].values
            iceThickness = ma.masked_values(iceThickness, 0)
            lon = mpasClimatology.remappedDataSet['lon'].values
            lat = mpasClimatology.remappedDataSet['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            for hemisphere in ['NH', 'SH']:

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

                obsClimatology = \
                    ObservationClimatology(
                        task=self,
                        fieldName=obsFieldName,
                        monthNames=months,
                        obsGridDescriptor=obsDescriptor,
                        comparisonGrid='latlon',
                        mappingFileSection='seaIceObservations',
                        mappingFileOption='seaIceClimatologyMappingFile',
                        mappingFilePrefix='map_obs_seaIce',
                        method=config.get('seaIceObservations',
                                          'interpolationMethod'))

                if obsClimatology.remappedDataSet is None:
                    obsClimatology.dataSet = xr.open_dataset(obsFileName)
                    obsClimatology.remap_and_write()

                obsIceThickness = \
                    obsClimatology.remappedDataSet['HI'].values

                # Mask thickness fields
                obsIceThickness = ma.masked_values(obsIceThickness, 0)
                if hemisphere == 'NH':
                    # Obs thickness should be nan above 86 (ICESat data)
                    obsIceThickness[latTarg > 86] = ma.masked
                    plotProjection = 'npstere'
                else:
                    plotProjection = 'spstere'

                difference = iceThickness - obsIceThickness

                title = '{} ({}, years {:04d}-{:04d})'.format(
                    subtitle, months, mpasClimatology.startYear,
                    mpasClimatology.endYear)
                fileout = '{}/icethick{}_{}_{}_years{:04d}-{:04d}.png'.format(
                    plotsDirectory, hemisphere, mainRunName, months,
                    mpasClimatology.startYear, mpasClimatology.endYear)
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

    def _open_sea_ice_part(self, inputFileNames, startDate, endDate):  # {{{
        """
        Open part of a sea-ice data set between the given start and end date,
        used to cache a climatology of the data set.

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

        variableList = [self.fieldName]

        ds = open_multifile_dataset(
            fileNames=inputFileNames,
            calendar=self.calendar,
            config=self.config,
            simulationStartTime=self.simulationStartTime,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=self.variableMap,
            startDate=startDate,
            endDate=endDate)

        return ds  # }}}

    # }}}
# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python