# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
#

import xarray as xr
import numpy as np
import netCDF4
import os
from mpas_tools.ocean.moc import add_moc_southern_boundary_transects
from mpas_tools.io import write_netcdf

from mpas_analysis.shared.constants.constants import m3ps_to_Sv
from mpas_analysis.shared.plot import plot_vertical_section_comparison, \
    timeseries_analysis_plot, savefig

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_files_year_month, get_region_mask

from mpas_analysis.shared.io import open_mpas_dataset\

from mpas_analysis.shared.timekeeping.utility import days_to_datetime

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.climatology.climatology import \
    get_climatology_op_directory

from mpas_analysis.shared.regions import ComputeRegionMasksSubtask


class StreamfunctionMOC(AnalysisTask):
    """
    Computation and plotting of model meridional overturning circulation.
    Will eventually support:

        * MOC streamfunction, post-processed
        * MOC streamfunction, from MOC analysis member
        * MOC time series (max value at 24.5N), post-processed
        * MOC time series (max value at 24.5N), from MOC analysis member
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  mpas_tools.config.MpasConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
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

        maskSubtask = ComputeMOCMasksSubtask(self)
        self.add_subtask(maskSubtask)

        computeClimSubtask = ComputeMOCClimatologySubtask(
            self, mpasClimatologyTask, maskSubtask)
        plotClimSubtask = PlotMOCClimatologySubtask(self, controlConfig)
        plotClimSubtask.run_after(computeClimSubtask)

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.get('timeSeries', 'endYear')
        if endYear == 'end':
            # a valid end year wasn't found, so likely the run was not found,
            # perhaps because we're just listing analysis tasks
            endYear = startYear
        else:
            endYear = int(endYear)

        years = range(startYear, endYear + 1)

        # in the end, we'll combine all the time series into one, but we create
        # this task first so it's easier to tell it to run after all the
        # compute tasks
        combineTimeSeriesSubtask = CombineMOCTimeSeriesSubtask(
            self, startYears=years, endYears=years)

        # run one subtask per year
        for year in years:
            computeTimeSeriesSubtask = ComputeMOCTimeSeriesSubtask(
                self, startYear=year, endYear=year, maskSubtask=maskSubtask)
            combineTimeSeriesSubtask.run_after(computeTimeSeriesSubtask)

        plotTimeSeriesSubtask = PlotMOCTimeSeriesSubtask(self, controlConfig)
        plotTimeSeriesSubtask.run_after(combineTimeSeriesSubtask)


class ComputeMOCMasksSubtask(ComputeRegionMasksSubtask):
    """
    An analysis subtasks for computing cell masks and southern transects for
    MOC regions
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask):

        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = parentTask.config
        meshName = config.get('input', 'mpasMeshName')
        regionGroup = 'MOC Basins'

        subprocessCount = config.getint('execute', 'parallelTaskCount')

        # call the constructor from the base class (ComputeRegionMasksSubtask)
        super().__init__(
            parentTask, regionGroup=regionGroup, meshName=meshName,
            subprocessCount=subprocessCount,
            useMpasMaskCreator=False)

        self.maskAndTransectFileName = None

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError :
            If a restart file is not available from which to read mesh
            information or if no history files are available from which to
            compute the climatology in the desired time range.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the parent class
        super().setup_and_check()

        self.maskAndTransectFileName = get_region_mask(
            self.config, '{}_mocBasinsAndTransects{}.nc'.format(
                self.meshName, self.date))

    def run_task(self):
        """
        Compute the requested climatologies
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if os.path.exists(self.maskAndTransectFileName):
            return

        # call ComputeRegionMasksSubtask.run_task() first
        super().run_task()

        config = self.config

        dsMesh = xr.open_dataset(self.obsFileName)
        dsMask = xr.open_dataset(self.maskFileName)

        dsMasksAndTransects = add_moc_southern_boundary_transects(
            dsMask, dsMesh, logger=self.logger)

        write_netcdf(dsMasksAndTransects, self.maskAndTransectFileName,
                     char_dim_name='StrLen')
# }}}


class ComputeMOCClimatologySubtask(AnalysisTask):
    """
    Computation of a climatology of the model meridional overturning
    circulation.

    Attributes
    ----------

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, mpasClimatologyTask, maskSubtask):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        maskSubtask : mpas_analysis.ocean.streamfunction_moc.ComputeMOCMasksSubtask
            The subtask for computing MOC region masks that runs before this
            subtask
        """
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
        self.maskSubtask = maskSubtask
        self.run_after(maskSubtask)

        parentTask.add_subtask(self)

        self.includeBolus = None
        self.includeSubmesoscale = None

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        """
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

        self.usePostprocessing = config.getexpression(
            self.sectionName, 'usePostprocessingScript')

        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            variableList = \
                ['timeMonthly_avg_mocStreamvalLatAndDepth',
                 'timeMonthly_avg_mocStreamvalLatAndDepthRegion']
        else:
            variableList = ['timeMonthly_avg_normalVelocity',
                            'timeMonthly_avg_vertVelocityTop',
                            'timeMonthly_avg_layerThickness']

            # Add the bolus velocity if GM is enabled
            try:
                # the new name
                self.includeBolus = self.namelist.getbool('config_use_gm')
            except KeyError:
                # the old name
                self.includeBolus = self.namelist.getbool(
                    'config_use_standardgm')
            try:
                self.includeSubmesoscale = \
                    self.namelist.getbool('config_submesoscale_enable')
            except KeyError:
                # an old run without submesoscale
                self.includeSubmesoscale = False

            if self.includeBolus:
                variableList.extend(
                    ['timeMonthly_avg_normalGMBolusVelocity',
                     'timeMonthly_avg_vertGMBolusVelocityTop'])

            if self.includeSubmesoscale:
                variableList.extend(
                    ['timeMonthly_avg_normalMLEvelocity',
                     'timeMonthly_avg_vertMLEBolusVelocityTop'])


        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=['ANN'])

    def run_task(self):
        """
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not.
        """
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

    def _compute_moc_climo_analysismember(self):
        """compute mean MOC streamfunction from analysis member"""

        config = self.config

        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        outputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear,
            self.endYear)

        if os.path.exists(outputFileName):
            return

        regionNames = config.getexpression(self.sectionName, 'regionNames')
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

        dsMask = xr.open_dataset(self.maskSubtask.maskAndTransectFileName)
        regionIndices = {}
        for iRegion in range(dsMask.sizes['nRegions']):
            regionInFile = str(dsMask.regionNames[iRegion].values.astype('U'))
            region = regionInFile.replace('_MOC', '')
            regionIndices[region] = iRegion

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
                indRegion = regionIndices[region]
                mocVar = annualClimatology.avgMocStreamfunRegional
                mocTop = mocVar.isel(nRegions=indRegion).values
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

    def _compute_moc_climo_postprocess(self):
        """compute mean MOC streamfunction as a post-process"""

        config = self.config
        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        outputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, self.startYear,
            self.endYear)

        if os.path.exists(outputFileName):
            return

        dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
            refTopDepth, refLayerThickness, cellsOnEdge = \
            _load_mesh(self.runStreams)

        regionNames = config.getexpression(self.sectionName, 'regionNames')

        # Load basin region related variables and save them to dictionary
        mpasMeshName = config.get('input', 'mpasMeshName')

        masksFileName = self.maskSubtask.maskAndTransectFileName
        dictRegion = _build_region_mask_dict(
            masksFileName, regionNames, mpasMeshName, self.logger)

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

        # rename some variables for convenience
        annualClimatology = annualClimatology.rename(
            {'timeMonthly_avg_normalVelocity': 'avgNormalVelocity',
             'timeMonthly_avg_vertVelocityTop': 'avgVertVelocityTop',
             'timeMonthly_avg_layerThickness': 'layerThickness'})

        if self.includeBolus:
            annualClimatology['avgNormalVelocity'] = \
                annualClimatology['avgNormalVelocity'] + \
                annualClimatology['timeMonthly_avg_normalGMBolusVelocity']

            annualClimatology['avgVertVelocityTop'] = \
                annualClimatology['avgVertVelocityTop'] + \
                annualClimatology['timeMonthly_avg_vertGMBolusVelocityTop']

        if self.includeSubmesoscale:
            annualClimatology['avgNormalVelocity'] = \
                annualClimatology['avgNormalVelocity'] + \
                annualClimatology['timeMonthly_avg_normalMLEvelocity']

            annualClimatology['avgVertVelocityTop'] = \
                annualClimatology['avgVertVelocityTop'] + \
                annualClimatology['timeMonthly_avg_vertMLEBolusVelocityTop']

        # Convert to numpy arrays
        # (can result in a memory error for large array size)
        horizontalVel = annualClimatology.avgNormalVelocity.values
        verticalVel = annualClimatology.avgVertVelocityTop.values
        velArea = verticalVel * areaCell[:, np.newaxis]
        layerThickness = annualClimatology.layerThickness.values

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
                                                horizontalVel,
                                                layerThickness,
                                                cellsOnEdge)

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


class PlotMOCClimatologySubtask(AnalysisTask):
    """
    Computation of a climatology of the model meridional overturning
    circulation.
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, controlConfig):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
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

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        """
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

        self.regionNames = ['Global'] + config.getexpression(self.sectionName,
                                                             'regionNames')

        for region in self.regionNames:
            filePrefix = 'moc{}_{}_years{:04d}-{:04d}'.format(
                region, mainRunName,
                self.startYear, self.endYear)

            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[region] = filePrefix

    def run_task(self):
        """
        Plot the MOC climatology
        """
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
            outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            x = lat[region]
            z = depth
            regionMOC = moc[region]
            # Subset lat range
            minLat = config.getexpression('streamfunctionMOC{}'.format(region),
                                          'latBinMin')
            maxLat = config.getexpression('streamfunctionMOC{}'.format(region),
                                          'latBinMax')
            indLat = np.logical_and(x >= minLat, x <= maxLat)
            x = x.where(indLat, drop=True)
            regionMOC = regionMOC.where(indLat, drop=True)
            if self.controlConfig is None:
                refRegionMOC = None
                diff = None
            else:
                # the coords of the ref MOC won't necessarily match this MOC
                # so we need to interpolate
                refRegionMOC = _interp_moc(x, z, regionMOC, refLat[region],
                                           refDepth, refMOC[region])

                diff = regionMOC - refRegionMOC

            plot_vertical_section_comparison(
                config, regionMOC, refRegionMOC, diff, xCoords=x, zCoord=z,
                colorMapSectionName='streamfunctionMOC{}'.format(region),
                colorbarLabel=colorbarLabel,
                title=title,
                modelTitle=mainRunName,
                refTitle=refTitle,
                diffTitle=diffTitle,
                xlabels=xLabel,
                ylabel=yLabel,
                movingAveragePoints=movingAveragePointsClimatological,
                maxTitleLength=70)

            savefig(outFileName, config)

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
                imageCaption=caption)

    def _load_moc(self, config):
        """compute mean MOC streamfunction from analysis member"""

        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        inputFileName = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
            outputDirectory, startYear,
            endYear)

        # Read from file
        ds = xr.open_dataset(inputFileName)
        depth = ds['depth']
        lat = {}
        moc = {}
        for region in self.regionNames:
            lat[region] = ds['lat{}'.format(region)]
            moc[region] = ds['moc{}'.format(region)]
        return depth, lat, moc


class ComputeMOCTimeSeriesSubtask(AnalysisTask):
    """
    Computation of a time series of max Atlantic MOC at 26.5N.
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, startYear, endYear, maskSubtask):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : mpas_analysis.ocean.streamfunction_moc.StreamfunctionMOC
            The main task of which this is a subtask

        startYear : int
            The start year of the time series

        endYear : int
            The end year of the time series

        maskSubtask : mpas_analysis.ocean.streamfunction_moc.ComputeMOCMasksSubtask
            The subtask for computing MOC region masks that runs before this
            subtask
        """
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

        self.maskSubtask = maskSubtask
        self.run_after(maskSubtask)

        parentTask.add_subtask(self)
        self.startYear = startYear
        self.endYear = endYear

        self.includeBolus = None
        self.includeSubmesoscale = None

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        """
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

        self.usePostprocessing = config.getexpression(
            self.sectionName, 'usePostprocessingScript')

        if not self.usePostprocessing and self.mocAnalysisMemberEnabled:
            self.variableList = \
                ['timeMonthly_avg_mocStreamvalLatAndDepth',
                 'timeMonthly_avg_mocStreamvalLatAndDepthRegion']
        else:
            self.variableList = ['timeMonthly_avg_normalVelocity',
                                 'timeMonthly_avg_vertVelocityTop',
                                 'timeMonthly_avg_layerThickness']

            # Add the bolus velocity if GM is enabled
            try:
                # the new name
                self.includeBolus = self.namelist.getbool('config_use_gm')
            except KeyError:
                # the old name
                self.includeBolus = self.namelist.getbool(
                    'config_use_standardgm')

            try:
                self.includeSubmesoscale = \
                    self.namelist.getbool('config_submesoscale_enable')
            except KeyError:
                # an old run without submesoscale
                self.includeSubmesoscale = False

            if self.includeBolus:
                self.variableList.extend(
                    ['timeMonthly_avg_normalGMBolusVelocity',
                     'timeMonthly_avg_vertGMBolusVelocityTop'])

            if self.includeSubmesoscale:
                self.variableList.extend(
                    ['timeMonthly_avg_normalMLEvelocity',
                     'timeMonthly_avg_vertMLEBolusVelocityTop'])

    def run_task(self):
        """
        Process MOC analysis member data if available, or compute MOC at
        post-processing if not.
        """
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

    def _compute_moc_time_series_analysismember(self):
        """compute MOC time series from analysis member"""

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
        moc = None
        refTopDepth = None
        times = np.zeros(len(inputFiles))
        computed = np.zeros(len(inputFiles), bool)

        continueOutput = os.path.exists(outputFileName)
        if continueOutput:
            self.logger.info('   Read in previously computed MOC time series')
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=['mocAtlantic26',
                                                 'mocAtlantic'],
                                   startDate=self.startDate,
                                   endDate=self.endDate) as dsMOCIn:

                dsMOCIn.load()

                if moc is None:
                    sizes = dsMOCIn.sizes
                    moc = np.zeros((len(inputFiles), sizes['depth'],
                                    sizes['lat']))
                    refTopDepth = dsMOCIn.depth.values

                # first, copy all computed data
                for inIndex in range(dsMOCIn.dims['Time']):

                    mask = np.logical_and(
                        dsMOCIn.year[inIndex].values == years,
                        dsMOCIn.month[inIndex].values == months)

                    outIndex = np.where(mask)[0][0]

                    mocRegion[outIndex] = dsMOCIn.mocAtlantic26[inIndex]
                    moc[outIndex, :, :] = dsMOCIn.mocAtlantic[inIndex, :, :]
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

            if moc is None:
                sizes = dsLocal.sizes
                moc = np.zeros((len(inputFiles), sizes['nVertLevels']+1,
                                len(binBoundaryMocStreamfunction)))
                try:
                    restartFile = self.runStreams.readpath('restart')[0]
                except ValueError:
                    raise IOError('No MPAS-O restart file found: need at '
                                  'least one restart file for MOC calculation')
                with xr.open_dataset(restartFile) as dsRestart:
                    refBottomDepth = dsRestart.refBottomDepth.values
                nVertLevels = len(refBottomDepth)
                refTopDepth = np.zeros(nVertLevels + 1)
                refTopDepth[1:nVertLevels + 1] = refBottomDepth[0:nVertLevels]

            moc[timeIndex, 0:-1, :] = mocTop

        description = 'Max MOC Atlantic streamfunction nearest to RAPID ' \
            'Array latitude (26.5N)'

        descriptionAtl = 'Atlantic MOC streamfunction'

        dictionary = {
            'dims': ['Time', 'depth', 'lat'],
            'coords': {
                'Time': {
                    'dims': ('Time',),
                    'data': times,
                    'attrs': {'units': 'days since 0001-01-01'}},
                'year': {
                    'dims': ('Time',),
                    'data': years,
                    'attrs': {'units': 'year'}},
                'month': {
                    'dims': ('Time',),
                    'data': months,
                    'attrs': {'units': 'month'}},
                'lat': {
                    'dims': ('lat',),
                    'data': binBoundaryMocStreamfunction,
                    'attrs': {'units': 'degrees north'}},
                'depth': {
                    'dims': ('depth',),
                    'data': refTopDepth,
                    'attrs': {'units': 'meters'}}},
            'data_vars': {
                'mocAtlantic26': {
                    'dims': ('Time',),
                    'data': mocRegion,
                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                              'description': description}},
                'mocAtlantic': {
                    'dims': ('Time', 'depth', 'lat'),
                    'data': moc,
                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                              'description': descriptionAtl}}}}
        dsMOCTimeSeries = xr.Dataset.from_dict(dictionary)
        write_netcdf(dsMOCTimeSeries, outputFileName)

    def _compute_moc_time_series_postprocess(self):
        """compute MOC time series as a post-process"""

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
            refTopDepth, refLayerThickness, cellsOnEdge = \
            _load_mesh(self.runStreams)

        mpasMeshName = config.get('input', 'mpasMeshName')

        masksFileName = self.maskSubtask.maskAndTransectFileName
        dictRegion = _build_region_mask_dict(
            masksFileName, ['Atlantic'], mpasMeshName, self.logger)
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
        moc = np.zeros((len(inputFiles), nVertLevels+1, len(latBins)))
        times = np.zeros(len(inputFiles))
        computed = np.zeros(len(inputFiles), bool)

        continueOutput = os.path.exists(outputFileName)
        if continueOutput:
            self.logger.info('   Read in previously computed MOC time series')
            with open_mpas_dataset(fileName=outputFileName,
                                   calendar=self.calendar,
                                   timeVariableNames=None,
                                   variableList=['mocAtlantic26',
                                                 'mocAtlantic'],
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
                    moc[outIndex, :, :] = dsMOCIn.mocAtlantic[inIndex, :, :]
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

            # rename some variables for convenience
            dsLocal = dsLocal.rename(
                {'timeMonthly_avg_normalVelocity': 'avgNormalVelocity',
                 'timeMonthly_avg_vertVelocityTop': 'avgVertVelocityTop',
                 'timeMonthly_avg_layerThickness': 'layerThickness'})

            if self.includeBolus:
                dsLocal['avgNormalVelocity'] = \
                    dsLocal['avgNormalVelocity'] + \
                    dsLocal['timeMonthly_avg_normalGMBolusVelocity']

                dsLocal['avgVertVelocityTop'] = \
                    dsLocal['avgVertVelocityTop'] + \
                    dsLocal['timeMonthly_avg_vertGMBolusVelocityTop']

            if self.includeSubmesoscale:
                dsLocal['avgNormalVelocity'] = \
                    dsLocal['avgNormalVelocity'] + \
                    dsLocal['timeMonthly_avg_normalMLEvelocity']

                dsLocal['avgVertVelocityTop'] = \
                    dsLocal['avgVertVelocityTop'] + \
                    dsLocal['timeMonthly_avg_vertMLEBolusVelocityTop']

            horizontalVel = dsLocal.avgNormalVelocity.values
            verticalVel = dsLocal.avgVertVelocityTop.values
            velArea = verticalVel * areaCell[:, np.newaxis]
            layerThickness = dsLocal.layerThickness.values

            transportZ = _compute_transport(maxEdgesInTransect,
                                            transectEdgeGlobalIDs,
                                            transectEdgeMaskSigns,
                                            nVertLevels, dvEdge,
                                            horizontalVel,
                                            layerThickness,
                                            cellsOnEdge)
            mocTop = _compute_moc(latAtlantic, nVertLevels, latCell,
                                  regionCellMask, transportZ, velArea)
            moc[timeIndex, :, :] = mocTop
            mocRegion[timeIndex] = np.amax(mocTop[:, indlat26])

        description = 'Max MOC Atlantic streamfunction nearest to RAPID ' \
            'Array latitude (26.5N)'

        descriptionAtl = 'Atlantic MOC streamfunction'

        dictionary = {
            'dims': ['Time', 'depth', 'lat'],
            'coords': {
                'Time': {
                    'dims': ('Time',),
                    'data': times,
                    'attrs': {'units': 'days since 0001-01-01'}},
                'year': {
                    'dims': ('Time',),
                    'data': years,
                    'attrs': {'units': 'year'}},
                'month': {
                    'dims': ('Time',),
                    'data': months,
                    'attrs': {'units': 'month'}},
                'lat': {
                    'dims': ('lat',),
                    'data': latAtlantic,
                    'attrs': {'units': 'degrees north'}},
                'depth': {
                    'dims': ('depth',),
                    'data': refTopDepth,
                    'attrs': {'units': 'meters'}}},
            'data_vars': {
                'mocAtlantic26': {
                    'dims': ('Time',),
                    'data': mocRegion,
                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                              'description': description}},
                'mocAtlantic': {
                    'dims': ('Time', 'depth', 'lat'),
                    'data': moc,
                    'attrs': {'units': 'Sv (10^6 m^3/s)',
                              'description': descriptionAtl}}}}
        dsMOCTimeSeries = xr.Dataset.from_dict(dictionary)
        write_netcdf(dsMOCTimeSeries, outputFileName)


class CombineMOCTimeSeriesSubtask(AnalysisTask):
    """
    Combine individual time series of max Atlantic MOC at 26.5N into a single
    data set
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
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

    def run_task(self):
        """
        Plot the MOC time series
        """
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
                               combine='nested', decode_times=False)

        ds.load()

        write_netcdf(ds, outputFileName)


class PlotMOCTimeSeriesSubtask(AnalysisTask):
    """
    Plots a time series of max Atlantic MOC at 26.5N.
    """
    # Authors
    # -------
    # Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, parentTask, controlConfig):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : ``StreamfunctionMOC``
            The main task of which this is a subtask

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
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

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError
            if timeSeriesStatsMonthly is not enabled in the MPAS run
        """
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

    def run_task(self):
        """
        Plot the MOC time series
        """
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
        title = '{}\n{}'.format(r'Max Atlantic MOC at $26.5\degree$N',
                                mainRunName)
        filePrefix = self.filePrefix

        outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

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
        lineColors = [config.get('timeSeries', 'mainColor')]

        lineWidths = [2]
        legendText = [mainRunName]

        if self.controlConfig is not None:

            dsRefMOC = self._load_moc(self.controlConfig)
            fields.append(dsRefMOC.mocAtlantic26)
            lineColors.append(config.get('timeSeries', 'controlColor'))
            lineWidths.append(2)
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            legendText.append(controlRunName)

        timeseries_analysis_plot(config, fields, calendar=self.calendar,
                                 title=title, xlabel=xLabel, ylabel=yLabel,
                                 movingAveragePoints=movingAveragePoints,
                                 lineColors=lineColors, lineWidths=lineWidths,
                                 legendText=legendText,
                                 firstYearXTicks=firstYearXTicks,
                                 yearStrideXTicks=yearStrideXTicks,
                                 maxTitleLength=90)

        savefig(outFileName, config)

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

    def _load_moc(self, config):
        """compute mean MOC streamfunction from analysis member"""

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        inputFileName = '{}/moc/mocTimeSeries_{:04d}-{:04d}.nc'.format(
            outputDirectory, startYear, endYear)

        dsMOCTimeSeries = xr.open_dataset(inputFileName, decode_times=False)
        return dsMOCTimeSeries


def _load_mesh(runStreams):
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
    cellsOnEdge = ncFile.variables['cellsOnEdge'][:] - 1
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
        refTopDepth, refLayerThickness, cellsOnEdge


def _build_region_mask_dict(regionMaskFile, regionNames, mpasMeshName, logger):
    if not os.path.exists(regionMaskFile):
        raise IOError('Regional masking file {} for MOC calculation '
                      'does not exist'.format(regionMaskFile))

    dsMask = xr.open_dataset(regionMaskFile)
    dsMask.load()

    regionIndices = {}
    for iRegion in range(dsMask.sizes['nRegions']):
        regionInFile = str(dsMask.regionNames[iRegion].values.astype('U'))
        region = regionInFile.replace('_MOC', '')
        regionIndices[region] = iRegion

    dictRegion = {}
    for region in regionNames:
        logger.info('\n  Reading region and transect mask for '
                    '{}...'.format(region))
        iRegion = regionIndices[region]
        maxEdgesInTransect = dsMask.sizes['maxEdgesInTransect']
        transectEdgeMaskSigns = \
            dsMask.transectEdgeMaskSigns.isel(nTransects=iRegion).values
        transectEdgeGlobalIDs = \
            dsMask.transectEdgeGlobalIDs.isel(nTransects=iRegion).values
        regionCellMask = \
            dsMask.regionCellMasks.isel(nRegions=iRegion).values

        indRegion = np.where(regionCellMask == 1)
        dictRegion[region] = {
            'indices': indRegion,
            'cellMask': regionCellMask,
            'maxEdgesInTransect': maxEdgesInTransect,
            'transectEdgeMaskSigns': transectEdgeMaskSigns,
            'transectEdgeGlobalIDs': transectEdgeGlobalIDs}

    return dictRegion


def _compute_transport(maxEdgesInTransect, transectEdgeGlobalIDs,
                       transectEdgeMaskSigns, nz, dvEdge,
                       horizontalVel, layerThickness, cellsOnEdge):
    """compute mass transport across southern transect of ocean basin"""

    transportZEdge = np.zeros([nz, maxEdgesInTransect])
    for i in range(maxEdgesInTransect):
        if transectEdgeGlobalIDs[i] == 0:
            break
        # subtract 1 because of python 0-indexing
        iEdge = transectEdgeGlobalIDs[i] - 1
        coe0 = cellsOnEdge[iEdge, 0]
        coe1 = cellsOnEdge[iEdge, 1]
        layerThicknessEdge = 0.5*(layerThickness[coe0, :] +
                                  layerThickness[coe1, :])
        transportZEdge[:, i] = horizontalVel[iEdge, :] * \
            transectEdgeMaskSigns[iEdge, np.newaxis] * \
            dvEdge[iEdge, np.newaxis] * \
            layerThicknessEdge[np.newaxis, :]
    transportZ = transportZEdge.sum(axis=1)
    return transportZ


def _compute_moc(latBins, nz, latCell, regionCellMask, transportZ,
                 velArea):
    """compute meridionally integrated MOC streamfunction"""

    mocTop = np.zeros([np.size(latBins), nz + 1])
    mocSouthBottomUp = - transportZ[::-1].cumsum()
    mocTop[0, 0:nz] = mocSouthBottomUp[::-1]
    for iLat in range(1, np.size(latBins)):
        indlat = np.logical_and(np.logical_and(
            regionCellMask == 1, latCell >= latBins[iLat - 1]),
            latCell < latBins[iLat])
        mocTop[iLat, :] = mocTop[iLat - 1, :] + \
            velArea[indlat, :].sum(axis=0)
    # convert m^3/s to Sverdrup
    mocTop = mocTop * m3ps_to_Sv
    mocTop = mocTop.T
    return mocTop


def _interp_moc(x, z, regionMOC, refX, refZ, refMOC):
    x = x.values
    z = z.values
    dims = regionMOC.dims
    regionMOC = regionMOC.values
    refX = refX.values
    refZ = refZ.values
    refMOC = refMOC.values

    nz, nx = regionMOC.shape
    refNz, refNx = refMOC.shape
    temp = np.zeros((refNz, nx))
    for zIndex in range(refNz):
        temp[zIndex, :] = np.interp(
            x, refX, refMOC[zIndex, :],
            left=np.nan, right=np.nan)
    refRegionMOC = np.zeros((nz, nx))
    for xIndex in range(nx):
        refRegionMOC[:, xIndex] = np.interp(
            z, refZ, temp[:, xIndex],
            left=np.nan, right=np.nan)

    return xr.DataArray(dims=dims, data=refRegionMOC)
