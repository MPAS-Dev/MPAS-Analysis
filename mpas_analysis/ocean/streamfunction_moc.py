# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import numpy as np
import netCDF4
import os

from mpas_analysis.shared.constants.constants import m3ps_to_Sv
from mpas_analysis.shared.plot import plot_vertical_section_comparison, \
    timeseries_analysis_plot

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_files_year_month

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf

from mpas_analysis.shared.timekeeping.utility import days_to_datetime

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.climatology.climatology import \
    get_climatology_op_directory


class StreamfunctionMOC(AnalysisTask):  # {{{
    '''
    Computation and plotting of model meridional overturning circulation.
    Will eventually support:

        * MOC streamfunction, post-processed
        * MOC streamfunction, from MOC analysis member
        * MOC time series (max value at 24.5N), post-processed
        * MOC time series (max value at 24.5N), from MOC analysis member
    '''
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(StreamfunctionMOC, self).__init__(
            config=config,
            taskName='streamfunctionMOC',
            componentName='ocean',
            tags=['streamfunction', 'moc', 'climatology', 'timeSeries',
                  'publicObs'])

        computeClimSubtask = ComputeMOCClimatologySubtask(
            self, mpasClimatologyTask)
        plotClimSubtask = PlotMOCClimatologySubtask(self, controlConfig)
        plotClimSubtask.run_after(computeClimSubtask)

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        years = range(startYear, endYear + 1)

        # in the end, we'll combine all the time series into one, but we create
        # this task first so it's easier to tell it to run after all the
        # compute tasks
        combineTimeSeriesSubtask = CombineMOCTimeSeriesSubtask(
            self, startYears=years, endYears=years)

        # run one subtask per year
        for year in years:
            computeTimeSeriesSubtask = ComputeMOCTimeSeriesSubtask(
                self, startYear=year, endYear=year)
            combineTimeSeriesSubtask.run_after(computeTimeSeriesSubtask)

        plotTimeSeriesSubtask = PlotMOCTimeSeriesSubtask(self, controlConfig)
        plotTimeSeriesSubtask.run_after(combineTimeSeriesSubtask)

        # }}}
    # }}}


class ComputeMOCClimatologySubtask(AnalysisTask):  # {{{
    '''
    Computation of a climatology of the model meridional overturning
    circulation.

    Attributes
    ----------

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    '''
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, mpasClimatologyTask):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeMOCClimatologySubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeMOCClimatology')

        self.mpasClimatologyTask = mpasClimatologyTask
        self.run_after(mpasClimatologyTask)

        parentTask.add_subtask(self)
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ComputeMOCClimatologySubtask, self).setup_and_check()

        self.startYear = self.mpasClimatologyTask.startYear
        self.startDate = self.mpasClimatologyTask.startDate
        self.endYear = self.mpasClimatologyTask.endYear
        self.endDate = self.mpasClimatologyTask.endDate

        config = self.config

        self.mocAnalysisMemberEnabled = self.check_analysis_enabled(
            analysisOptionName='config_am_mocstreamfunction_enable',
            raiseException=False)

        self.sectionName = 'streamfunctionMOC'

        self.usePostprocessing = config.getExpression(
            self.sectionName, 'usePostprocessingScript')

        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            variableList = \
                ['timeMonthly_avg_mocStreamvalLatAndDepth',
                 'timeMonthly_avg_mocStreamvalLatAndDepthRegion']
        else:
            variableList = ['timeMonthly_avg_normalVelocity',
                            'timeMonthly_avg_vertVelocityTop']

            # Add the bolus velocity if GM is enabled
            self.includeBolus = self.namelist.getbool('config_use_standardgm')
            if self.includeBolus:
                variableList.extend(
                    ['timeMonthly_avg_normalGMBolusVelocity',
                     'timeMonthly_avg_vertGMBolusVelocityTop'])

        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=['ANN'])

        # }}}

    def run_task(self):  # {{{
        '''
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not.
        '''
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("Computing climatology of Meridional Overturning "
                         "Circulation (MOC)...")

        # **** Compute MOC ****
        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            self._compute_moc_climo_analysismember()
        else:
            self._compute_moc_climo_postprocess()

        # }}}

    def _compute_moc_climo_analysismember(self):  # {{{
        '''compute mean MOC streamfunction from analysis member'''

        config = self.config

        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        outputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear,
            self.endYear)

        if os.path.exists(outputFileName):
            return

        regionNames = config.getExpression(self.sectionName, 'regionNames')
        regionNames.append('Global')

        # Read in depth and bin latitudes
        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least '
                          'one for MHT calcuation')

        with xr.open_dataset(restartFileName) as dsRestart:
            refBottomDepth = dsRestart.refBottomDepth.values

        nVertLevels = len(refBottomDepth)
        refLayerThickness = np.zeros(nVertLevels)
        refLayerThickness[0] = refBottomDepth[0]
        refLayerThickness[1:nVertLevels] = \
            refBottomDepth[1:nVertLevels] - refBottomDepth[0:nVertLevels - 1]

        refZMid = refBottomDepth - 0.5 * refLayerThickness

        binBoundaryMocStreamfunction = None
        # first try timeSeriesStatsMonthly for bin boundaries, then try
        # mocStreamfunctionOutput stream as a backup option
        for streamName in ['timeSeriesStatsMonthlyOutput',
                           'mocStreamfunctionOutput']:
            try:
                inputFileName = self.historyStreams.readpath(streamName)[0]
            except ValueError:
                raise IOError('At least one file from stream {} is needed '
                              'to compute MOC'.format(streamName))

            with xr.open_dataset(inputFileName) as ds:
                if 'binBoundaryMocStreamfunction' in ds.data_vars:
                    binBoundaryMocStreamfunction = \
                        ds.binBoundaryMocStreamfunction.values
                    break

        if binBoundaryMocStreamfunction is None:
            raise ValueError('Could not find binBoundaryMocStreamfunction in '
                             'either timeSeriesStatsMonthlyOutput or '
                             'mocStreamfunctionOutput streams')

        binBoundaryMocStreamfunction = np.rad2deg(binBoundaryMocStreamfunction)

        # Compute and plot annual climatology of MOC streamfunction
        self.logger.info('\n  Compute climatology of MOC streamfunction...')
        self.logger.info('   Load data...')

        climatologyFileName = self.mpasClimatologyTask.get_file_name(
            season='ANN')
        annualClimatology = xr.open_dataset(climatologyFileName)
        if 'Time' in annualClimatology.dims:
            annualClimatology = annualClimatology.isel(Time=0)

        # rename some variables for convenience
        annualClimatology = annualClimatology.rename(
            {'timeMonthly_avg_mocStreamvalLatAndDepth':
                'avgMocStreamfunGlobal',
             'timeMonthly_avg_mocStreamvalLatAndDepthRegion':
                 'avgMocStreamfunRegional'})

        # Create dictionary for MOC climatology (NB: need this form
        # in order to convert it to xarray dataset later in the script)
        depth = refZMid
        lat = {}
        moc = {}
        for region in regionNames:
            self.logger.info('   Compute {} MOC...'.format(region))
            if region == 'Global':
                mocTop = annualClimatology.avgMocStreamfunGlobal.values
            else:
                # hard-wire region=0 (Atlantic) for now
                indRegion = 0
                mocVar = annualClimatology.avgMocStreamfunRegional
                mocTop = mocVar[indRegion, :, :].values
            # Store computed MOC to dictionary
            lat[region] = binBoundaryMocStreamfunction
            moc[region] = mocTop

        # Save to file
        self.logger.info('   Save global and regional MOC to file...')
        ncFile = netCDF4.Dataset(outputFileName, mode='w')
        # create dimensions
        ncFile.createDimension('nz', nVertLevels)
        for region in regionNames:
            latBins = lat[region]
            mocTop = moc[region]
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
        depthVar = ncFile.createVariable('depth', 'f4', ('nz',))
        depthVar.description = 'depth'
        depthVar.units = 'meters'
        depthVar[:] = depth
        ncFile.close()
        # }}}

    def _compute_moc_climo_postprocess(self):  # {{{
        '''compute mean MOC streamfunction as a post-process'''

        config = self.config
        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        outputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear,
            self.endYear)

        if os.path.exists(outputFileName):
            return

        dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness = _load_mesh(self.runStreams)

        regionNames = config.getExpression(self.sectionName, 'regionNames')

        # Load basin region related variables and save them to dictionary
        mpasMeshName = config.get('input', 'mpasMeshName')
        regionMaskDirectory = build_config_full_path(config, 'diagnostics',
                                                     'regionMasksubDirectory')

        dictRegion = _build_region_mask_dict(regionNames, regionMaskDirectory,
                                             mpasMeshName, self.logger)

        # Add Global regionCellMask=1 everywhere to make the algorithm
        # for the global moc similar to that of the regional moc
        dictRegion['Global'] = {
            'cellMask': np.ones(np.size(latCell))}
        regionNames.append('Global')

        # Compute and plot annual climatology of MOC streamfunction
        self.logger.info('\n  Compute post-processed climatological of MOC '
                         'streamfunction...')

        self.logger.info('   Load data...')

        climatologyFileName = self.mpasClimatologyTask.get_file_name(
            season='ANN')
        annualClimatology = xr.open_dataset(climatologyFileName)
        if 'Time' in annualClimatology.dims:
            annualClimatology = annualClimatology.isel(Time=0)

        if self.includeBolus:
            annualClimatology['avgNormalVelocity'] = \
                annualClimatology['timeMonthly_avg_normalVelocity'] + \
                annualClimatology['timeMonthly_avg_normalGMBolusVelocity']

            annualClimatology['avgVertVelocityTop'] = \
                annualClimatology['timeMonthly_avg_vertVelocityTop'] + \
                annualClimatology['timeMonthly_avg_vertGMBolusVelocityTop']
        else:
            # rename some variables for convenience
            annualClimatology = annualClimatology.rename(
                {'timeMonthly_avg_normalVelocity': 'avgNormalVelocity',
                 'timeMonthly_avg_vertVelocityTop': 'avgVertVelocityTop'})

        # Convert to numpy arrays
        # (can result in a memory error for large array size)
        horizontalVel = annualClimatology.avgNormalVelocity.values
        verticalVel = annualClimatology.avgVertVelocityTop.values
        velArea = verticalVel * areaCell[:, np.newaxis]

        # Create dictionary for MOC climatology (NB: need this form
        # in order to convert it to xarray dataset later in the script)
        depth = refTopDepth
        lat = {}
        moc = {}
        for region in regionNames:
            self.logger.info('   Compute {} MOC...'.format(region))
            self.logger.info('    Compute transport through region '
                             'southern transect...')
            if region == 'Global':
                transportZ = np.zeros(nVertLevels)
            else:
                maxEdgesInTransect = \
                    dictRegion[region]['maxEdgesInTransect']
                transectEdgeGlobalIDs = \
                    dictRegion[region]['transectEdgeGlobalIDs']
                transectEdgeMaskSigns = \
                    dictRegion[region]['transectEdgeMaskSigns']
                transportZ = _compute_transport(maxEdgesInTransect,
                                                transectEdgeGlobalIDs,
                                                transectEdgeMaskSigns,
                                                nVertLevels, dvEdge,
                                                refLayerThickness,
                                                horizontalVel)

            regionCellMask = dictRegion[region]['cellMask']
            latBinSize = \
                config.getfloat('streamfunctionMOC{}'.format(region),
                                'latBinSize')
            if region == 'Global':
                latBins = np.arange(-90.0, 90.1, latBinSize)
            else:
                indRegion = dictRegion[region]['indices']
                latBins = latCell[indRegion]
                latBins = np.arange(np.amin(latBins),
                                    np.amax(latBins) + latBinSize,
                                    latBinSize)
            mocTop = _compute_moc(latBins, nVertLevels, latCell,
                                  regionCellMask, transportZ, velArea)

            # Store computed MOC to dictionary
            lat[region] = latBins
            moc[region] = mocTop

        # Save to file
        self.logger.info('   Save global and regional MOC to file...')
        ncFile = netCDF4.Dataset(outputFileName, mode='w')
        # create dimensions
        ncFile.createDimension('nz', len(refTopDepth))
        for region in regionNames:
            latBins = lat[region]
            mocTop = moc[region]
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
        depthVar = ncFile.createVariable('depth', 'f4', ('nz',))
        depthVar.description = 'depth'
        depthVar.units = 'meters'
        depthVar[:] = depth
        ncFile.close()
        # }}}

    # }}}


class PlotMOCClimatologySubtask(AnalysisTask):  # {{{
    '''
    Computation of a climatology of the model meridional overturning
    circulation.
    '''
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, controlConfig):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotMOCClimatologySubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plotMOCClimatology')

        parentTask.add_subtask(self)

        self.controlConfig = controlConfig
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotMOCClimatologySubtask, self).setup_and_check()

        config = self.config

        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')

        self.sectionName = 'streamfunctionMOC'

        self.xmlFileNames = []
        self.filePrefixes = {}

        mainRunName = config.get('runs', 'mainRunName')

        self.regionNames = ['Global'] + config.getExpression(self.sectionName,
                                                             'regionNames')

        for region in self.regionNames:
            filePrefix = 'moc{}_{}_years{:04d}-{:04d}'.format(
                region, mainRunName,
                self.startYear, self.endYear)

            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[region] = filePrefix

        # }}}

    def run_task(self):  # {{{
        '''
        Plot the MOC climatology
        '''
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("\nPlotting streamfunction of Meridional Overturning "
                         "Circulation (MOC)...")

        config = self.config

        depth, lat, moc = self._load_moc(config)

        if self.controlConfig is None:
            refTitle = None
            diffTitle = None
        else:
            refDepth, refLat, refMOC = self._load_moc(self.controlConfig)
            refTitle = self.controlConfig.get('runs', 'mainRunName')
            diffTitle = 'Main - Control'

        # **** Plot MOC ****
        # Define plotting variables
        mainRunName = config.get('runs', 'mainRunName')
        movingAveragePointsClimatological = config.getint(
            self.sectionName, 'movingAveragePointsClimatological')
        colorbarLabel = '[Sv]'
        xLabel = 'latitude [deg]'
        yLabel = 'depth [m]'

        for region in self.regionNames:
            self.logger.info('   Plot climatological {} MOC...'.format(region))
            title = '{} MOC (ANN, years {:04d}-{:04d})'.format(
                region, self.startYear,
                self.endYear)
            filePrefix = self.filePrefixes[region]
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            x = lat[region]
            z = depth
            regionMOC = moc[region]
            # Subset lat range
            minLat = config.getExpression('streamfunctionMOC{}'.format(region),
                                          'latBinMin')
            maxLat = config.getExpression('streamfunctionMOC{}'.format(region),
                                          'latBinMax')
            indLat = np.logical_and(x >= minLat, x <= maxLat)
            x = x[indLat]
            regionMOC = regionMOC[:, indLat]
            if self.controlConfig is None:
                refRegionMOC = None
                diff = None
            else:
                # the coords of the ref MOC won't necessarily match this MOC
                # so we need to interpolate
                nz, nx = regionMOC.shape
                refNz, refNx = refMOC[region].shape
                temp = np.zeros((refNz, nx))
                for zIndex in range(refNz):
                    temp[zIndex, :] = np.interp(
                        x, refLat[region], refMOC[region][zIndex, :],
                        left=np.nan, right=np.nan)
                refRegionMOC = np.zeros((nz, nx))
                for xIndex in range(nx):
                    refRegionMOC[:, xIndex] = np.interp(
                        depth, refDepth, temp[:, xIndex],
                        left=np.nan, right=np.nan)

                diff = regionMOC - refRegionMOC

            plot_vertical_section_comparison(
                config, x, z, regionMOC, refRegionMOC, diff,
                fileout=figureName,
                colorMapSectionName='streamfunctionMOC{}'.format(region),
                cbarLabel=colorbarLabel,
                title=title,
                modelTitle=mainRunName,
                refTitle=refTitle,
                diffTitle=diffTitle,
                xlabel=xLabel,
                ylabel=yLabel,
                N=movingAveragePointsClimatological)

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

    def _load_moc(self, config):  # {{{
        '''compute mean MOC streamfunction from analysis member'''

        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        inputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, startYear,
            endYear)

        # Read from file
        ncFile = netCDF4.Dataset(inputFileName, mode='r')
        depth = ncFile.variables['depth'][:]
        lat = {}
        moc = {}
        for region in self.regionNames:
            lat[region] = ncFile.variables['lat{}'.format(region)][:]
            moc[region] = \
                ncFile.variables['moc{}'.format(region)][:, :]
        ncFile.close()
        return depth, lat, moc  # }}}

    # }}}


class ComputeMOCTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Computation of a time series of max Atlantic MOC at 26.5N.
    '''
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, startYear, endYear):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeMOCTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='computeMOCTimeSeries_{:04d}-{:04d}'.format(
                startYear, endYear))

        parentTask.add_subtask(self)
        self.startYear = startYear
        self.endYear = endYear
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ComputeMOCTimeSeriesSubtask, self).setup_and_check()

        config = self.config

        self.mocAnalysisMemberEnabled = self.check_analysis_enabled(
            analysisOptionName='config_am_mocstreamfunction_enable',
            raiseException=False)

        self.sectionName = 'streamfunctionMOC'

        self.usePostprocessing = config.getExpression(
            self.sectionName, 'usePostprocessingScript')

        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            self.variableList = \
                ['timeMonthly_avg_mocStreamvalLatAndDepth',
                 'timeMonthly_avg_mocStreamvalLatAndDepthRegion']
        else:
            self.variableList = ['timeMonthly_avg_normalVelocity',
                                 'timeMonthly_avg_vertVelocityTop']

            # Add the bolus velocity if GM is enabled
            self.includeBolus = self.namelist.getbool('config_use_standardgm')
            if self.includeBolus:
                self.variableList.extend(
                    ['timeMonthly_avg_normalGMBolusVelocity',
                     'timeMonthly_avg_vertGMBolusVelocityTop'])
        # }}}

    def run_task(self):  # {{{
        '''
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not.
        '''
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("\nCompute time series of Meridional Overturning "
                         "Circulation (MOC)...")

        self.startDate = '{:04d}-01-01_00:00:00'.format(self.startYear)
        self.endDate = '{:04d}-12-31_23:59:59'.format(self.endYear)

        # **** Compute MOC ****
        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            self._compute_moc_time_series_analysismember()
        else:
            self._compute_moc_time_series_postprocess()
        # }}}

    def _compute_moc_time_series_analysismember(self):  # {{{
        '''compute MOC time series from analysis member'''

        # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
        self.logger.info('\n  Compute Atlantic MOC time series from analysis '
                         'member...')
        self.logger.info('   Load data...')

        outputDirectory = '{}/moc/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'))
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileName = '{}/mocTimeSeries_{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear, self.endYear)

        streamName = 'timeSeriesStatsMonthlyOutput'

        # Get bin latitudes and index of 26.5N
        binBoundaryMocStreamfunction = None
        # first try timeSeriesStatsMonthly for bin boundaries, then try
        # mocStreamfunctionOutput stream as a backup option
        for streamName in ['timeSeriesStatsMonthlyOutput',
                           'mocStreamfunctionOutput']:
            try:
                inputFileName = self.historyStreams.readpath(streamName)[0]
            except ValueError:
                raise IOError('At least one file from stream {} is needed '
                              'to compute MOC'.format(streamName))

            with xr.open_dataset(inputFileName) as ds:
                if 'binBoundaryMocStreamfunction' in ds.data_vars:
                    binBoundaryMocStreamfunction = \
                        ds.binBoundaryMocStreamfunction.values
                    break

        if binBoundaryMocStreamfunction is None:
            raise ValueError('Could not find binBoundaryMocStreamfunction in '
                             'either timeSeriesStatsMonthlyOutput or '
                             'mocStreamfunctionOutput streams')

        binBoundaryMocStreamfunction = np.rad2deg(binBoundaryMocStreamfunction)
        dLat = binBoundaryMocStreamfunction - 26.5
        indlat26 = np.where(np.abs(dLat) == np.amin(np.abs(dLat)))

        inputFiles = sorted(self.historyStreams.readpath(
            streamName, startDate=self.startDate,
            endDate=self.endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        mocRegion = np.zeros(len(inputFiles))
        times = np.zeros(len(inputFiles))
        computed = np.zeros(len(inputFiles), bool)

        continueOutput = os.path.exists(outputFileName)
        if continueOutput:
            self.logger.info('   Read in previously computed MOC time series')
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=['mocAtlantic26'],
                                   startDate=self.startDate,
                                   endDate=self.endDate) as dsMOCIn:

                dsMOCIn.load()

                # first, copy all computed data
                for inIndex in range(dsMOCIn.dims['Time']):

                    mask = np.logical_and(
                        dsMOCIn.year[inIndex].values == years,
                        dsMOCIn.month[inIndex].values == months)

                    outIndex = np.where(mask)[0][0]

                    mocRegion[outIndex] = dsMOCIn.mocAtlantic26[inIndex]
                    times[outIndex] = dsMOCIn.Time[inIndex]
                    computed[outIndex] = True

                if np.all(computed):
                    # no need to waste time writing out the data set again
                    return dsMOCIn

        for timeIndex, fileName in enumerate(inputFiles):
            if computed[timeIndex]:
                continue

            dsLocal = open_mpas_dataset(
                fileName=fileName,
                calendar=self.calendar,
                variableList=self.variableList,
                startDate=self.startDate,
                endDate=self.endDate)
            dsLocal = dsLocal.isel(Time=0)
            time = dsLocal.Time.values
            times[timeIndex] = time
            date = days_to_datetime(time, calendar=self.calendar)

            self.logger.info('     date: {:04d}-{:02d}'.format(date.year,
                                                               date.month))

            # hard-wire region=0 (Atlantic) for now
            indRegion = 0
            mocVar = dsLocal.timeMonthly_avg_mocStreamvalLatAndDepthRegion
            mocTop = mocVar[indRegion, :, :].values
            mocRegion[timeIndex] = np.amax(mocTop[:, indlat26])

        description = 'Max MOC Atlantic streamfunction nearest to RAPID ' \
            'Array latitude (26.5N)'

        dictonary = {'dims': ['Time'],
                     'coords': {'Time':
                                {'dims': ('Time'),
                                 'data': times,
                                 'attrs': {'units': 'days since 0001-01-01'}},
                                'year':
                                {'dims': ('Time'),
                                 'data': years,
                                 'attrs': {'units': 'year'}},
                                'month':
                                {'dims': ('Time'),
                                 'data': months,
                                 'attrs': {'units': 'month'}}},
                     'data_vars': {'mocAtlantic26':
                                   {'dims': ('Time'),
                                    'data': mocRegion,
                                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                                              'description': description}}}}
        dsMOCTimeSeries = xr.Dataset.from_dict(dictonary)
        write_netcdf(dsMOCTimeSeries, outputFileName)
        # }}}

    def _compute_moc_time_series_postprocess(self):  # {{{
        '''compute MOC time series as a post-process'''

        config = self.config

        # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
        self.logger.info('\n  Compute and/or plot post-processed Atlantic MOC '
                         'time series...')
        self.logger.info('   Load data...')

        outputDirectory = '{}/moc/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'))
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileName = '{}/mocTimeSeries_{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear, self.endYear)

        dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness = _load_mesh(self.runStreams)

        mpasMeshName = config.get('input', 'mpasMeshName')
        regionMaskDirectory = build_config_full_path(config, 'diagnostics',
                                                     'regionMasksubDirectory')
        dictRegion = _build_region_mask_dict(['Atlantic'], regionMaskDirectory,
                                             mpasMeshName, self.logger)
        dictRegion = dictRegion['Atlantic']

        latBinSize = config.getfloat('streamfunctionMOCAtlantic',
                                     'latBinSize')
        indRegion = dictRegion['indices']
        latBins = latCell[indRegion]
        latBins = np.arange(np.amin(latBins),
                            np.amax(latBins) + latBinSize,
                            latBinSize)
        latAtlantic = latBins
        dLat = latAtlantic - 26.5
        indlat26 = np.where(np.abs(dLat) == np.amin(np.abs(dLat)))

        maxEdgesInTransect = dictRegion['maxEdgesInTransect']
        transectEdgeGlobalIDs = dictRegion['transectEdgeGlobalIDs']
        transectEdgeMaskSigns = dictRegion['transectEdgeMaskSigns']
        regionCellMask = dictRegion['cellMask']

        streamName = 'timeSeriesStatsMonthlyOutput'
        inputFiles = sorted(self.historyStreams.readpath(
            streamName, startDate=self.startDate,
            endDate=self.endDate, calendar=self.calendar))

        years, months = get_files_year_month(inputFiles,
                                             self.historyStreams,
                                             'timeSeriesStatsMonthlyOutput')

        mocRegion = np.zeros(len(inputFiles))
        times = np.zeros(len(inputFiles))
        computed = np.zeros(len(inputFiles), bool)

        continueOutput = os.path.exists(outputFileName)
        if continueOutput:
            self.logger.info('   Read in previously computed MOC time series')
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=['mocAtlantic26'],
                                   startDate=self.startDate,
                                   endDate=self.endDate) as dsMOCIn:

                dsMOCIn.load()

                # first, copy all computed data
                for inIndex in range(dsMOCIn.dims['Time']):

                    mask = np.logical_and(
                        dsMOCIn.year[inIndex].values == years,
                        dsMOCIn.month[inIndex].values == months)

                    outIndex = np.where(mask)[0][0]

                    mocRegion[outIndex] = dsMOCIn.mocAtlantic26[inIndex]
                    times[outIndex] = dsMOCIn.Time[inIndex]
                    computed[outIndex] = True

                if np.all(computed):
                    # no need to waste time writing out the data set again
                    return dsMOCIn

        for timeIndex, fileName in enumerate(inputFiles):
            if computed[timeIndex]:
                continue

            dsLocal = open_mpas_dataset(
                fileName=fileName,
                calendar=self.calendar,
                variableList=self.variableList,
                startDate=self.startDate,
                endDate=self.endDate)
            dsLocal = dsLocal.isel(Time=0)
            time = dsLocal.Time.values
            times[timeIndex] = time
            date = days_to_datetime(time, calendar=self.calendar)

            self.logger.info('     date: {:04d}-{:02d}'.format(date.year,
                                                               date.month))

            if self.includeBolus:
                dsLocal['avgNormalVelocity'] = \
                    dsLocal['timeMonthly_avg_normalVelocity'] + \
                    dsLocal['timeMonthly_avg_normalGMBolusVelocity']

                dsLocal['avgVertVelocityTop'] = \
                    dsLocal['timeMonthly_avg_vertVelocityTop'] + \
                    dsLocal['timeMonthly_avg_vertGMBolusVelocityTop']
            else:
                # rename some variables for convenience
                dsLocal = dsLocal.rename(
                    {'timeMonthly_avg_normalVelocity': 'avgNormalVelocity',
                     'timeMonthly_avg_vertVelocityTop': 'avgVertVelocityTop'})

            horizontalVel = dsLocal.avgNormalVelocity.values
            verticalVel = dsLocal.avgVertVelocityTop.values
            velArea = verticalVel * areaCell[:, np.newaxis]
            transportZ = _compute_transport(maxEdgesInTransect,
                                            transectEdgeGlobalIDs,
                                            transectEdgeMaskSigns,
                                            nVertLevels, dvEdge,
                                            refLayerThickness,
                                            horizontalVel)
            mocTop = _compute_moc(latAtlantic, nVertLevels, latCell,
                                  regionCellMask, transportZ, velArea)
            mocRegion[timeIndex] = np.amax(mocTop[:, indlat26])

        description = 'Max MOC Atlantic streamfunction nearest to RAPID ' \
            'Array latitude (26.5N)'

        dictonary = {'dims': ['Time'],
                     'coords': {'Time':
                                {'dims': ('Time'),
                                 'data': times,
                                 'attrs': {'units': 'days since 0001-01-01'}},
                                'year':
                                {'dims': ('Time'),
                                 'data': years,
                                 'attrs': {'units': 'year'}},
                                'month':
                                {'dims': ('Time'),
                                 'data': months,
                                 'attrs': {'units': 'month'}}},
                     'data_vars': {'mocAtlantic26':
                                   {'dims': ('Time'),
                                    'data': mocRegion,
                                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                                              'description': description}}}}
        dsMOCTimeSeries = xr.Dataset.from_dict(dictonary)
        write_netcdf(dsMOCTimeSeries, outputFileName)
        # }}}
    # }}}


class CombineMOCTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Combine individual time series of max Atlantic MOC at 26.5N into a single
    data set
    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(CombineMOCTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='combineMOCTimeSeries')

        parentTask.add_subtask(self)
        self.startYears = startYears
        self.endYears = endYears
        # }}}

    def run_task(self):  # {{{
        '''
        Plot the MOC time series
        '''
        # Authors
        # -------
        # Xylar Asay-Davis
        outputDirectory = '{}/moc/'.format(
            build_config_full_path(self.config, 'output',
                                   'timeseriesSubdirectory'))
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outputFileNames = []
        for startYear, endYear in zip(self.startYears, self.endYears):
            outputFileName = '{}/mocTimeSeries_{:04d}-{:04d}.nc'.format(
                outputDirectory, startYear, endYear)
            outputFileNames.append(outputFileName)

        outputFileName = '{}/mocTimeSeries_{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYears[0], self.endYears[-1])

        if outputFileName in outputFileNames:
            # don't try to write to read from and write to the same file
            return

        ds = xr.open_mfdataset(outputFileNames, concat_dim='Time',
                               decode_times=False)

        write_netcdf(ds, outputFileName)  # }}}
    # }}}


class PlotMOCTimeSeriesSubtask(AnalysisTask):  # {{{
    '''
    Plots a time series of max Atlantic MOC at 26.5N.
    '''
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, controlConfig):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotMOCTimeSeriesSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plotMOCTimeSeries')

        parentTask.add_subtask(self)

        self.controlConfig = controlConfig
        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(PlotMOCTimeSeriesSubtask, self).setup_and_check()

        config = self.config

        self.sectionName = 'streamfunctionMOC'

        mainRunName = config.get('runs', 'mainRunName')

        filePrefix = 'mocTimeseries_{}'.format(mainRunName)
        self.xmlFileNames = ['{}/{}.xml'.format(self.plotsDirectory,
                                                filePrefix)]
        self.filePrefix = filePrefix

        # }}}

    def run_task(self):  # {{{
        '''
        Plot the MOC time series
        '''
        # Authors
        # -------
        # Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

        self.logger.info("\nPlotting time series of Meridional Overturning "
                         "Circulation (MOC)...")

        config = self.config

        dsMOCTimeSeries = self._load_moc(config)

        # **** Plot MOC ****
        # Define plotting variables
        mainRunName = config.get('runs', 'mainRunName')
        movingAveragePoints = config.getint(self.sectionName,
                                            'movingAveragePoints')

        # Plot time series
        self.logger.info('   Plot time series of max Atlantic MOC at 26.5N...')
        xLabel = 'Time [years]'
        yLabel = '[Sv]'
        title = r'Max Atlantic MOC at $26.5\degree$N\n {}'.format(mainRunName)
        filePrefix = self.filePrefix

        figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

        if config.has_option(self.taskName, 'firstYearXTicks'):
            firstYearXTicks = config.getint(self.taskName,
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option(self.taskName, 'yearStrideXTicks'):
            yearStrideXTicks = config.getint(self.taskName,
                                             'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        fields = [dsMOCTimeSeries.mocAtlantic26]
        lineColors = ['k']
        lineWidths = [2]
        legendText = [mainRunName]

        if self.controlConfig is not None:

            dsRefMOC = self._load_moc(self.controlConfig)
            fields.append(dsRefMOC.mocAtlantic26)
            lineColors.append('r')
            lineWidths.append(2)
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            legendText.append(controlRunName)

        timeseries_analysis_plot(config, fields,
                                 movingAveragePoints, title,
                                 xLabel, yLabel, figureName,
                                 calendar=self.calendar, lineColors=lineColors,
                                 lineWidths=lineWidths,
                                 legendText=legendText,
                                 firstYearXTicks=firstYearXTicks,
                                 yearStrideXTicks=yearStrideXTicks)

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
            imageCaption=caption)

        # }}}

    def _load_moc(self, config):  # {{{
        '''compute mean MOC streamfunction from analysis member'''

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        inputFileName = '{}/moc/mocTimeSeries_{:04d}-{:04d}.nc'.format(
            outputDirectory, startYear, endYear)

        dsMOCTimeSeries = xr.open_dataset(inputFileName, decode_times=False)
        return dsMOCTimeSeries  # }}}

    # }}}


def _load_mesh(runStreams):  # {{{
    # Load mesh related variables
    try:
        restartFile = runStreams.readpath('restart')[0]
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
    refTopDepth = np.zeros(nVertLevels + 1)
    refTopDepth[1:nVertLevels + 1] = refBottomDepth[0:nVertLevels]
    refLayerThickness = np.zeros(nVertLevels)
    refLayerThickness[0] = refBottomDepth[0]
    refLayerThickness[1:nVertLevels] = \
        (refBottomDepth[1:nVertLevels] -
         refBottomDepth[0:nVertLevels - 1])

    return dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        refTopDepth, refLayerThickness
    # }}}


def _build_region_mask_dict(regionNames, regionMaskDirectory, mpasMeshName,
                            logger):  # {{{

    regionMaskFile = '{}/{}_SingleRegionAtlanticWTransportTransects_' \
                     'masks.nc'.format(regionMaskDirectory, mpasMeshName)

    if not os.path.exists(regionMaskFile):
        raise IOError('Regional masking file {} for MOC calculation '
                      'does not exist'.format(regionMaskFile))
    iRegion = 0
    dictRegion = {}
    for region in regionNames:
        logger.info('\n  Reading region and transect mask for '
                    '{}...'.format(region))
        ncFileRegional = netCDF4.Dataset(regionMaskFile, mode='r')
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
        dictRegion[region] = {
            'indices': indRegion,
            'cellMask': regionCellMask,
            'maxEdgesInTransect': maxEdgesInTransect,
            'transectEdgeMaskSigns': transectEdgeMaskSigns,
            'transectEdgeGlobalIDs': transectEdgeGlobalIDs}

    return dictRegion  # }}}


def _compute_transport(maxEdgesInTransect, transectEdgeGlobalIDs,
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


def _compute_moc(latBins, nz, latCell, regionCellMask, transportZ,
                 velArea):  # {{{
    '''compute meridionally integrated MOC streamfunction'''

    mocTop = np.zeros([np.size(latBins), nz + 1])
    mocTop[0, range(1, nz + 1)] = transportZ.cumsum()
    for iLat in range(1, np.size(latBins)):
        indlat = np.logical_and(np.logical_and(
            regionCellMask == 1, latCell >= latBins[iLat - 1]),
            latCell < latBins[iLat])
        mocTop[iLat, :] = mocTop[iLat - 1, :] + \
            velArea[indlat, :].sum(axis=0)
    # convert m^3/s to Sverdrup
    mocTop = mocTop * m3ps_to_Sv
    mocTop = mocTop.T
    return mocTop  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
