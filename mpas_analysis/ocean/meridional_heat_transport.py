import xarray as xr
import numpy as np
import netCDF4
import os
import warnings

from ..shared.plot.plotting import plot_vertical_section,\
    setup_colormap, plot_1D

from ..shared.io.utility import build_config_full_path

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared import AnalysisTask
from ..shared.html import write_image_xml


class MeridionalHeatTransport(AnalysisTask):  # {{{
    '''
    Plot meridional heat transport from the analysis member output.

    Attributes
    ----------

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    Authors
    -------
    Mark Petersen, Milena Veneziani, Xylar Asay-Davis
    '''

    def __init__(self, config, mpasClimatologyTask):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        Authors
        -------
        Xylar Asay-Davis

        '''
        # first, call the constructor from the base class (AnalysisTask)
        super(MeridionalHeatTransport, self).__init__(
            config=config,
            taskName='meridionalHeatTransport',
            componentName='ocean',
            tags=['climatology'])

        self.mpasClimatologyTask = mpasClimatologyTask
        self.run_after(mpasClimatologyTask)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        ValueError: if myArg has an invalid value

        Authors
        -------
        Mark Petersen, Milena Veneziani, Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(MeridionalHeatTransport, self).setup_and_check()

        self.startYear = self.mpasClimatologyTask.startYear
        self.startDate = self.mpasClimatologyTask.startDate
        self.endYear = self.mpasClimatologyTask.endYear
        self.endDate = self.mpasClimatologyTask.endDate

        config = self.config

        self.check_analysis_enabled(
            analysisOptionName='config_am_meridionalheattransport_enable',
            raiseException=True)

        # Later, we will read in depth and MHT latitude points
        # from mpaso.hist.am.meridionalHeatTransport.*.nc
        mhtFiles = self.historyStreams.readpath(
                'meridionalHeatTransportOutput')
        if len(mhtFiles) == 0:
            raise IOError('No MPAS-O MHT history file found: need at least '
                          'one ')

        self.mhtFile = mhtFiles[0]

        self.simulationStartTime = get_simulation_start_time(self.runStreams)

        self.sectionName = 'meridionalHeatTransport'

        # Read in obs file information
        compareWithObs = config.getboolean(self.sectionName,
                                           'compareWithObservations')
        self.observationsFile = None
        if compareWithObs:
            observationsDirectory = build_config_full_path(
                config, 'oceanObservations', 'mhtSubdirectory')
            observationsFile = config.get(self.sectionName, 'observationData')
            observationsFile = '{}/{}'.format(observationsDirectory,
                                              observationsFile)
            if os.path.exists(observationsFile):
                self.observationsFile = observationsFile
            else:
                warnings.warn('No MHT observations file found: skip plotting '
                              'obs')

        mainRunName = self.config.get('runs', 'mainRunName')

        variableList = ['timeMonthly_avg_meridionalHeatTransportLat',
                        'timeMonthly_avg_meridionalHeatTransportLatZ']

        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=['ANN'])

        self.xmlFileNames = []
        self.filePrefixes = {}

        prefixes = ['mht']
        if config.getboolean(self.sectionName, 'plotVerticalSection'):
            prefixes.append('mhtZ')

        for prefix in prefixes:
            filePrefix = '{}_{}_years{:04d}-{:04d}'.format(
                    prefix, mainRunName,
                    self.startYear, self.endYear)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[prefix] = filePrefix

        # }}}

    def run_task(self):  # {{{
        """
        Process MHT analysis member data if available.
        Plots MHT as:
           1D function of latitude
           2D function of latitude and depth

        Authors
        -------
        Mark Petersen, Milena Veneziani, Xylar Asay-Davis
        """
        self.logger.info("\nPlotting meridional heat transport (MHT)...")

        config = self.config

        movingAveragePoints = config.getint('meridionalHeatTransport',
                                            'movingAveragePoints')

        # Read in depth and MHT latitude points
        # Latitude is from binBoundaryMerHeatTrans written in
        #  mpaso.hist.am.meridionalHeatTransport.*.nc
        # Depth is from refZMid, also in
        # mpaso.hist.am.meridionalHeatTransport.*.nc

        self.logger.info('  Read in depth and latitude...')
        ncFile = netCDF4.Dataset(self.mhtFile, mode='r')
        # reference depth [m]
        refZMid = ncFile.variables['refZMid'][:]
        refBottomDepth = ncFile.variables['refBottomDepth'][:]
        binBoundaryMerHeatTrans = \
            ncFile.variables['binBoundaryMerHeatTrans'][:]
        binBoundaryMerHeatTrans = np.rad2deg(binBoundaryMerHeatTrans)
        ncFile.close()

        nVertLevels = len(refBottomDepth)
        refLayerThickness = np.zeros(nVertLevels)
        refLayerThickness[0] = refBottomDepth[0]
        refLayerThickness[1:nVertLevels] = (refBottomDepth[1:nVertLevels] -
                                            refBottomDepth[0:nVertLevels-1])

        ######################################################################
        # Mark P Note: Currently only supports global MHT.
        # Need to add variables merHeatTransLatRegion and
        # merHeatTransLatZRegion
        # These are not computed by default in ACME right now.
        # Then we will need to add another section for regions with a loop
        # over number of regions.
        ######################################################################

        self.logger.info('\n   Plotting global meridional heat transport')

        self.logger.info('   Load data...')

        climatologyFileName = self.mpasClimatologyTask.get_file_name(
                season='ANN')

        annualClimatology = xr.open_dataset(climatologyFileName)
        annualClimatology = annualClimatology.isel(Time=0)

        # **** Plot MHT ****
        # Define plotting variables
        mainRunName = config.get('runs', 'mainRunName')
        xLimGlobal = config.getExpression(self.sectionName, 'xLimGlobal')
        depthLimGlobal = config.getExpression(self.sectionName,
                                              'depthLimGlobal')

        self.logger.info('   Plot global MHT...')
        # Plot 1D MHT (zonally averaged, depth integrated)
        x = binBoundaryMerHeatTrans
        y = annualClimatology.timeMonthly_avg_meridionalHeatTransportLat
        xLabel = 'latitude [deg]'
        yLabel = 'meridional heat transport [PW]'
        title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
                 self.startYear, self.endYear, mainRunName)
        filePrefix = self.filePrefixes['mht']
        figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        if self.observationsFile is not None:
            # Load in observations
            dsObs = xr.open_dataset(self.observationsFile)
            xObs = dsObs.LATITUDE
            ncepGlobal = dsObs.GLOBALNCEP_ADJUSTED
            ncepErrGlobal = dsObs.GLOBALNCEP_ERR
            ecmwfGlobal = dsObs.GLOBALECMWF_ADJUSTED
            ecmwfErrGlobal = dsObs.GLOBALECMWF_ERR

            lineColors = ['r', 'b', 'g']
            lineWidths = [1.6, 1.2, 1.2]
            legendText = ['model', 'NCEP', 'ECMWF']
            plot_1D(config, [x, xObs, xObs],
                    [y, ncepGlobal, ecmwfGlobal],
                    [None, ncepErrGlobal, ecmwfErrGlobal],
                    lineColors, lineWidths, legendText,
                    title, xLabel, yLabel, figureName,
                    xLim=xLimGlobal)
        else:
            lineColors = ['r']
            lineWidths = [1.6]
            legendText = [None]
            plot_1D(config, [x], [y], [None],
                    lineColors, lineWidths, legendText,
                    title, xLabel, yLabel, figureName,
                    xLim=xLimGlobal)

        self._write_xml(filePrefix)

        if config.getboolean(self.sectionName, 'plotVerticalSection'):
            # Plot 2D MHT (zonally integrated)

            # normalize 2D MHT by layer thickness
            MHTLatZVar = \
                annualClimatology.timeMonthly_avg_meridionalHeatTransportLatZ
            MHTLatZ = MHTLatZVar.values.T[:, :]
            for k in range(nVertLevels):
                MHTLatZ[k, :] = MHTLatZ[k, :]/refLayerThickness[k]

            x = binBoundaryMerHeatTrans
            y = refZMid
            z = MHTLatZ
            xLabel = 'latitude [deg]'
            yLabel = 'depth [m]'
            title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
                     self.startYear, self.endYear, mainRunName)
            filePrefix = self.filePrefixes['mhtZ']
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
            colorbarLabel = '[PW/m]'
            contourLevels = config.getExpression(self.sectionName,
                                                 'contourLevelsGlobal',
                                                 usenumpyfunc=True)
            (colormapName, colorbarLevels) = setup_colormap(config,
                                                            self.sectionName,
                                                            suffix='Global')
            plot_vertical_section(config, x, y, z,
                                  colormapName, colorbarLevels,
                                  contourLevels, colorbarLabel,
                                  title, xLabel, yLabel, figureName,
                                  xLim=xLimGlobal, yLim=depthLimGlobal,
                                  invertYAxis=False, N=movingAveragePoints)

            self._write_xml(filePrefix)

        # }}}

    def _write_xml(self, filePrefix):  # {{{
        caption = 'Meridional Heat Transport'
        write_image_xml(
            config=self.config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Meridional Heat Transport',
            groupLink='mht',
            imageDescription=caption,
            imageCaption=caption)  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
