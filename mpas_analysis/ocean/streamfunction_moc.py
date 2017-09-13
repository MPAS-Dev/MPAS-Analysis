# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import netCDF4
import os
from functools import partial

from ..shared.constants.constants import m3ps_to_Sv
from ..shared.plot.plotting import plot_vertical_section,\
    timeseries_analysis_plot, setup_colormap

from ..shared.io.utility import build_config_full_path, make_directories

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    days_to_datetime

from ..shared.climatology.climatology \
    import update_climatology_bounds_from_file_names, \
    compute_climatologies_with_ncclimo

from ..shared.analysis_task import AnalysisTask

from ..shared.time_series import cache_time_series
from ..shared.html import write_image_xml


class StreamfunctionMOC(AnalysisTask):  # {{{
    '''
    Computation and plotting of model meridional overturning circulation.
    Will eventually support:
      * MOC streamfunction, post-processed (currently supported)
      * MOC streamfunction, from MOC analysis member
      * MOC time series (max value at 24.5N), post-processed
      * MOC time series (max value at 24.5N), from MOC analysis member

    Authors
    -------
    Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis
    '''

    def __init__(self, config):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis

        '''
        # first, call the constructor from the base class (AnalysisTask)
        super(StreamfunctionMOC, self).__init__(
            config=config,
            taskName='streamfunctionMOC',
            componentName='ocean',
            tags=['streamfunction', 'moc', 'climatology', 'timeSeries'])

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(StreamfunctionMOC, self).setup_and_check()

        config = self.config

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        self.mocAnalysisMemberEnabled = self.check_analysis_enabled(
            analysisOptionName='config_am_mocstreamfunction_enable',
            raiseException=False)

        # Get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        #   First a list necessary for the streamfunctionMOC climatology
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.startDateClimo = config.get('climatology', 'startDate')
        self.endDateClimo = config.get('climatology', 'endDate')
        self.inputFilesClimo = \
            self.historyStreams.readpath(streamName,
                                         startDate=self.startDateClimo,
                                         endDate=self.endDateClimo,
                                         calendar=self.calendar)
        if len(self.inputFilesClimo) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDateClimo,
                                       self.endDateClimo))

        self.simulationStartTime = get_simulation_start_time(self.runStreams)

        changed, self.startYearClimo, self.endYearClimo, self.startDateClimo, \
            self.endDateClimo = update_climatology_bounds_from_file_names(
                    self.inputFilesClimo, self.config)

        #   Then a list necessary for the streamfunctionMOC Atlantic timeseries
        self.startDateTseries = config.get('timeSeries', 'startDate')
        self.endDateTseries = config.get('timeSeries', 'endDate')
        self.inputFilesTseries = \
            self.historyStreams.readpath(streamName,
                                         startDate=self.startDateTseries,
                                         endDate=self.endDateTseries,
                                         calendar=self.calendar)
        if len(self.inputFilesTseries) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDateTseries,
                                       self.endDateTseries))

        self.startYearTseries = config.getint('timeSeries', 'startYear')
        self.endYearTseries = config.getint('timeSeries', 'endYear')

        self.sectionName = 'streamfunctionMOC'

        self.xmlFileNames = []
        self.filePrefixes = {}

        mainRunName = config.get('runs', 'mainRunName')

        regions = ['Global'] + config.getExpression(self.sectionName,
                                                    'regionNames')

        for region in regions:
            filePrefix = 'moc{}_{}_years{:04d}-{:04d}'.format(
                    region, mainRunName,
                    self.startYearClimo, self.endYearClimo)

            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[region] = filePrefix

        filePrefix = 'mocTimeseries_{}'.format(mainRunName)
        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    filePrefix))
        self.filePrefixes['timeSeries'] = filePrefix

        # }}}

    def run(self):  # {{{
        '''
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not. Plots streamfunction climatolgoical sections
        as well as time series of max Atlantic MOC at 26.5N (latitude of
        RAPID MOC Array).

        Authors
        -------
        Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis
        '''

        print "\nPlotting streamfunction of Meridional Overturning " \
              "Circulation (MOC)..."

        print '\n  List of files for climatologies:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFilesClimo[0]),
                  os.path.basename(self.inputFilesClimo[-1]))

        print '\n  List of files for time series:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFilesTseries[0]),
                  os.path.basename(self.inputFilesTseries[-1]))

        config = self.config

        # **** Compute MOC ****
        # Check whether MOC Analysis Member is enabled
        if self.mocAnalysisMemberEnabled:
            # Add a moc_analisysMember_processing
            print '*** MOC Analysis Member is on ***'
            # (mocDictClimo, mocDictTseries) = \
            #     self._compute_moc_analysismember(config, streams, calendar,
            #                                      sectionName, dictClimo,
            #                                      dictTseries)
        else:
            self._compute_velocity_climatologies()
            self._compute_moc_climo_postprocess()
            dsMOCTimeSeries = self._compute_moc_time_series_postprocess()

        # **** Plot MOC ****
        # Define plotting variables
        mainRunName = config.get('runs', 'mainRunName')
        movingAveragePoints = config.getint(self.sectionName,
                                            'movingAveragePoints')
        colorbarLabel = '[Sv]'
        xLabel = 'latitude [deg]'
        yLabel = 'depth [m]'

        for region in self.regionNames:
            print '   Plot climatological {} MOC...'.format(region)
            title = '{} MOC (ANN, years {:04d}-{:04d})\n {}'.format(
                     region, self.startYearClimo,
                     self.endYearClimo,
                     mainRunName)
            filePrefix = self.filePrefixes[region]
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
            contourLevels = \
                config.getExpression(self.sectionName,
                                     'contourLevels{}'.format(region),
                                     usenumpyfunc=True)
            (colormapName, colorbarLevels) = setup_colormap(config,
                                                            self.sectionName,
                                                            suffix=region)

            x = self.lat[region]
            y = self.depth
            z = self.moc[region]
            plot_vertical_section(config, x, y, z, colormapName,
                                  colorbarLevels, contourLevels, colorbarLabel,
                                  title, xLabel, yLabel, figureName)

            caption = '{} Meridional Overturning Streamfunction'.format(region)
            write_image_xml(
                config=config,
                filePrefix=filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup='Meridional Overturning Streamfunction',
                groupLink='moc',
                thumbnailDescription=region,
                imageDescription=caption,
                imageCaption=caption)  # }}}

        # Plot time series
        print '   Plot time series of max Atlantic MOC at 26.5N...'
        xLabel = 'Time [years]'
        yLabel = '[Sv]'
        title = 'Max Atlantic MOC at $26.5^\circ$N\n {}'.format(mainRunName)
        filePrefix = self.filePrefixes['timeSeries']

        figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        timeseries_analysis_plot(config, [dsMOCTimeSeries.mocAtlantic26],
                                 movingAveragePoints, title,
                                 xLabel, yLabel, figureName,
                                 lineStyles=['k-'], lineWidths=[1.5],
                                 calendar=self.calendar)

        caption = u'Time Series of maximum Meridional Overturning ' \
                  u'Circulation at 26.5°N'
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Meridional Overturning Streamfunction',
            groupLink='moc',
            thumbnailDescription='Time Series',
            imageDescription=caption,
            imageCaption=caption)  # }}}

        # }}}

    def _load_mesh(self):  # {{{
        # Load mesh related variables
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for MOC calculation')
        ncFile = netCDF4.Dataset(restartFile, mode='r')
        dvEdge = ncFile.variables['dvEdge'][:]
        areaCell = ncFile.variables['areaCell'][:]
        refBottomDepth = ncFile.variables['refBottomDepth'][:]
        latCell = np.rad2deg(ncFile.variables['latCell'][:])
        ncFile.close()
        nVertLevels = len(refBottomDepth)
        refTopDepth = np.zeros(nVertLevels+1)
        refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]
        refLayerThickness = np.zeros(nVertLevels)
        refLayerThickness[0] = refBottomDepth[0]
        refLayerThickness[1:nVertLevels] = \
            (refBottomDepth[1:nVertLevels] -
             refBottomDepth[0:nVertLevels-1])

        return dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness
        # }}}

    def _compute_velocity_climatologies(self):  # {{{
        '''compute yearly velocity climatologies and cache them'''

        variableList = ['timeMonthly_avg_normalVelocity',
                        'timeMonthly_avg_vertVelocityTop']

        config = self.config

        outputRoot = build_config_full_path(config, 'output',
                                            'mpasClimatologySubdirectory')

        outputDirectory = '{}/meanVelocity'.format(outputRoot)

        self.velClimoFile = \
            '{}/mpaso_ANN_{:04d}01_{:04d}12_climo.nc'.format(
                    outputDirectory, self.startYearClimo,
                    self.endYearClimo)

        if not os.path.exists(self.velClimoFile):
            make_directories(outputDirectory)

            compute_climatologies_with_ncclimo(
                    config=config,
                    inDirectory=self.historyDirectory,
                    outDirectory=outputDirectory,
                    startYear=self.startYearClimo,
                    endYear=self.endYearClimo,
                    variableList=variableList,
                    modelName='mpaso',
                    decemberMode='sdd')
        # }}}

    def _compute_moc_climo_postprocess(self):  # {{{

        '''compute mean MOC streamfunction as a post-process'''

        config = self.config

        dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness = self._load_mesh()

        self.regionNames = config.getExpression(self.sectionName,
                                                'regionNames')

        # Load basin region related variables and save them to dictionary
        # NB: The following will need to change with new regional mapping files
        regionMaskFiles = config.get(self.sectionName, 'regionMaskFiles')
        if not os.path.exists(regionMaskFiles):
            raise IOError('Regional masking file for MOC calculation '
                          'does not exist')
        iRegion = 0
        self.dictRegion = {}
        for region in self.regionNames:
            print '\n  Reading region and transect mask for ' \
                '{}...'.format(region)
            ncFileRegional = netCDF4.Dataset(regionMaskFiles, mode='r')
            maxEdgesInTransect = \
                ncFileRegional.dimensions['maxEdgesInTransect'].size
            transectEdgeMaskSigns = \
                ncFileRegional.variables['transectEdgeMaskSigns'][:, iRegion]
            transectEdgeGlobalIDs = \
                ncFileRegional.variables['transectEdgeGlobalIDs'][iRegion, :]
            regionCellMask = \
                ncFileRegional.variables['regionCellMasks'][:, iRegion]
            ncFileRegional.close()
            iRegion += 1

            indRegion = np.where(regionCellMask == 1)
            self.dictRegion[region] = {
                'indices': indRegion,
                'cellMask': regionCellMask,
                'maxEdgesInTransect': maxEdgesInTransect,
                'transectEdgeMaskSigns': transectEdgeMaskSigns,
                'transectEdgeGlobalIDs': transectEdgeGlobalIDs}
        # Add Global regionCellMask=1 everywhere to make the algorithm
        # for the global moc similar to that of the regional moc

        self.dictRegion['Global'] = {
                'cellMask': np.ones(np.size(latCell))}
        self.regionNames.append('Global')

        # Compute and plot annual climatology of MOC streamfunction
        print '\n  Compute and/or plot post-processed MOC climatological '\
              'streamfunction...'
        outputDirectory = build_config_full_path(config, 'output',
                                                 'mpasClimatologySubdirectory')

        make_directories(outputDirectory)

        outputFileClimo = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
                           outputDirectory, self.startYearClimo,
                           self.endYearClimo)
        if not os.path.exists(outputFileClimo):
            print '   Load data...'

            annualClimatology = xr.open_dataset(self.velClimoFile)
            # rename some variables for convenience
            annualClimatology = annualClimatology.rename(
                    {'timeMonthly_avg_normalVelocity': 'avgNormalVelocity',
                     'timeMonthly_avg_vertVelocityTop': 'avgVertVelocityTop'})
            annualClimatology = annualClimatology.isel(Time=0)

            # Convert to numpy arrays
            # (can result in a memory error for large array size)
            horizontalVel = annualClimatology.avgNormalVelocity.values
            verticalVel = annualClimatology.avgVertVelocityTop.values
            velArea = verticalVel * areaCell[:, np.newaxis]

            # Create dictionary for MOC climatology (NB: need this form
            # in order to convert it to xarray dataset later in the script)
            self.depth = refTopDepth
            self.lat = {}
            self.moc = {}
            for region in self.regionNames:
                print '   Compute {} MOC...'.format(region)
                print '    Compute transport through region southern ' \
                    'transect...'
                if region == 'Global':
                    transportZ = np.zeros(nVertLevels)
                else:
                    maxEdgesInTransect = \
                        self.dictRegion[region]['maxEdgesInTransect']
                    transectEdgeGlobalIDs = \
                        self.dictRegion[region]['transectEdgeGlobalIDs']
                    transectEdgeMaskSigns = \
                        self.dictRegion[region]['transectEdgeMaskSigns']
                    transportZ = self._compute_transport(maxEdgesInTransect,
                                                         transectEdgeGlobalIDs,
                                                         transectEdgeMaskSigns,
                                                         nVertLevels, dvEdge,
                                                         refLayerThickness,
                                                         horizontalVel)

                regionCellMask = self.dictRegion[region]['cellMask']
                latBinSize = \
                    config.getExpression(self.sectionName,
                                         'latBinSize{}'.format(region))
                if region == 'Global':
                    latBins = np.arange(-90.0, 90.1, latBinSize)
                else:
                    indRegion = self.dictRegion[region]['indices']
                    latBins = latCell[indRegion]
                    latBins = np.arange(np.amin(latBins),
                                        np.amax(latBins)+latBinSize,
                                        latBinSize)
                mocTop = self._compute_moc(latBins, nVertLevels, latCell,
                                           regionCellMask, transportZ, velArea)

                # Store computed MOC to dictionary
                self.lat[region] = latBins
                self.moc[region] = mocTop

            # Save to file
            print '   Save global and regional MOC to file...'
            ncFile = netCDF4.Dataset(outputFileClimo, mode='w')
            # create dimensions
            ncFile.createDimension('nz', len(refTopDepth))
            for region in self.regionNames:
                latBins = self.lat[region]
                mocTop = self.moc[region]
                ncFile.createDimension('nx{}'.format(region), len(latBins))
                # create variables
                x = ncFile.createVariable('lat{}'.format(region), 'f4',
                                          ('nx{}'.format(region),))
                x.description = 'latitude bins for MOC {}'\
                                ' streamfunction'.format(region)
                x.units = 'degrees (-90 to 90)'
                y = ncFile.createVariable('moc{}'.format(region), 'f4',
                                          ('nz', 'nx{}'.format(region)))
                y.description = 'MOC {} streamfunction, annual'\
                                ' climatology'.format(region)
                y.units = 'Sv (10^6 m^3/s)'
                # save variables
                x[:] = latBins
                y[:, :] = mocTop
            depth = ncFile.createVariable('depth', 'f4', ('nz',))
            depth.description = 'depth'
            depth.units = 'meters'
            depth[:] = refTopDepth
            ncFile.close()
        else:
            # Read from file
            print '   Read previously computed MOC streamfunction from file...'
            ncFile = netCDF4.Dataset(outputFileClimo, mode='r')
            self.depth = ncFile.variables['depth'][:]
            self.lat = {}
            self.moc = {}
            for region in self.regionNames:
                self.lat[region] = ncFile.variables['lat{}'.format(region)][:]
                self.moc[region] = \
                    ncFile.variables['moc{}'.format(region)][:, :]
            ncFile.close()
        # }}}

    def _compute_moc_time_series_postprocess(self):  # {{{
        '''compute MOC time series as a post-process'''

        # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
        print '\n  Compute and/or plot post-processed Atlantic MOC '\
              'time series...'
        print '   Load data...'

        config = self.config

        self.simulationStartTime = get_simulation_start_time(self.runStreams)
        variableList = ['timeMonthly_avg_normalVelocity',
                        'timeMonthly_avg_vertVelocityTop']

        dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness = self._load_mesh()

        if config.has_option(self.sectionName, 'maxChunkSize'):
            chunking = config.getExpression(self.sectionName, 'maxChunkSize')
        else:
            chunking = None

        ds = open_multifile_dataset(
            fileNames=self.inputFilesTseries,
            calendar=self.calendar,
            config=config,
            simulationStartTime=self.simulationStartTime,
            timeVariableName=['xtime_startMonthly', 'xtime_endMonthly'],
            variableList=variableList,
            startDate=self.startDateTseries,
            endDate=self.endDateTseries,
            chunking=chunking)
        latAtlantic = self.lat['Atlantic']
        dLat = latAtlantic - 26.5
        indlat26 = np.where(dLat == np.amin(np.abs(dLat)))

        dictRegion = self.dictRegion['Atlantic']
        maxEdgesInTransect = dictRegion['maxEdgesInTransect']
        transectEdgeGlobalIDs = dictRegion['transectEdgeGlobalIDs']
        transectEdgeMaskSigns = dictRegion['transectEdgeMaskSigns']
        regionCellMask = dictRegion['cellMask']

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileTseries = '{}/mocTimeSeries.nc'.format(outputDirectory)

        continueOutput = os.path.exists(outputFileTseries)
        if continueOutput:
            print '   Read in previously computed MOC time series'

        # add all the other arguments to the function
        comp_moc_part = partial(self._compute_moc_time_series_part, ds,
                                areaCell, latCell, indlat26,
                                maxEdgesInTransect, transectEdgeGlobalIDs,
                                transectEdgeMaskSigns,  nVertLevels, dvEdge,
                                refLayerThickness, latAtlantic, regionCellMask)

        dsMOCTimeSeries = cache_time_series(
            ds.Time.values,  comp_moc_part, outputFileTseries,
            self.calendar, yearsPerCacheUpdate=1,  printProgress=False)

        return dsMOCTimeSeries  # }}}

    def _compute_moc_time_series_part(self, ds, areaCell, latCell, indlat26,
                                      maxEdgesInTransect,
                                      transectEdgeGlobalIDs,
                                      transectEdgeMaskSigns, nVertLevels,
                                      dvEdge, refLayerThickness, latAtlantic,
                                      regionCellMask, timeIndices, firstCall):
        # computes a subset of the MOC time series

        if firstCall:
            print '   Process and save time series'

        times = ds.Time[timeIndices].values
        mocRegion = np.zeros(timeIndices.shape)

        for localIndex, timeIndex in enumerate(timeIndices):
            time = times[localIndex]
            dsLocal = ds.isel(Time=timeIndex)
            date = days_to_datetime(time, calendar=self.calendar)

            print '     date: {:04d}-{:02d}'.format(date.year, date.month)

            horizontalVel = dsLocal.timeMonthly_avg_normalVelocity.values
            verticalVel = dsLocal.timeMonthly_avg_vertVelocityTop.values
            velArea = verticalVel * areaCell[:, np.newaxis]
            transportZ = self._compute_transport(maxEdgesInTransect,
                                                 transectEdgeGlobalIDs,
                                                 transectEdgeMaskSigns,
                                                 nVertLevels, dvEdge,
                                                 refLayerThickness,
                                                 horizontalVel)
            mocTop = self._compute_moc(latAtlantic, nVertLevels, latCell,
                                       regionCellMask, transportZ, velArea)
            mocRegion[localIndex] = np.amax(mocTop[:, indlat26])

        description = 'Max MOC Atlantic streamfunction nearest to RAPID ' \
            'Array latitude (26.5N)'

        dictonary = {'dims': ['Time'],
                     'coords': {'Time':
                                {'dims': ('Time'),
                                 'data': times,
                                 'attrs': {'units': 'days since 0001-01-01'}}},
                     'data_vars': {'mocAtlantic26':
                                   {'dims': ('Time'),
                                    'data': mocRegion,
                                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                                              'description': description}}}}
        dsMOC = xr.Dataset.from_dict(dictonary)
        return dsMOC

    # def _compute_moc_analysismember(self):
    #
    #     return

    def _compute_transport(self, maxEdgesInTransect, transectEdgeGlobalIDs,
                           transectEdgeMaskSigns, nz, dvEdge,
                           refLayerThickness, horizontalVel):  # {{{

        '''compute mass transport across southern transect of ocean basin'''

        transportZEdge = np.zeros([nz, maxEdgesInTransect])
        for i in range(maxEdgesInTransect):
            if transectEdgeGlobalIDs[i] == 0:
                break
            # subtract 1 because of python 0-indexing
            iEdge = transectEdgeGlobalIDs[i] - 1
            transportZEdge[:, i] = horizontalVel[iEdge, :] * \
                transectEdgeMaskSigns[iEdge, np.newaxis] * \
                dvEdge[iEdge, np.newaxis] * \
                refLayerThickness[np.newaxis, :]
        transportZ = transportZEdge.sum(axis=1)
        return transportZ  # }}}

    def _compute_moc(self, latBins, nz, latCell, regionCellMask, transportZ,
                     velArea):  # {{{

        '''compute meridionally integrated MOC streamfunction'''

        mocTop = np.zeros([np.size(latBins), nz+1])
        mocTop[0, range(1, nz+1)] = transportZ.cumsum()
        for iLat in range(1, np.size(latBins)):
            indlat = np.logical_and(np.logical_and(
                         regionCellMask == 1, latCell >= latBins[iLat-1]),
                         latCell < latBins[iLat])
            mocTop[iLat, :] = mocTop[iLat-1, :] + \
                velArea[indlat, :].sum(axis=0)
        # convert m^3/s to Sverdrup
        mocTop = mocTop * m3ps_to_Sv
        mocTop = mocTop.T
        return mocTop  # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
