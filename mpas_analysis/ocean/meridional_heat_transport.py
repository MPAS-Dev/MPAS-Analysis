import xarray as xr
import numpy as np
import netCDF4
import os
import warnings

from ..shared.plot.plotting import plot_vertical_section,\
    setup_colormap, plot_1D

from ..shared.io.utility import build_config_full_path

from ..shared.climatology import MpasClimatology

from ..shared.analysis_task import AnalysisTask


class MeridionalHeatTransport(AnalysisTask):  # {{{
    '''
    Plot meridional heat transport from the analysis member output.

    Authors
    -------
    Mark Petersen, Milena Veneziani, Xylar Asay-Davis
    '''

    @classmethod
    def create_tasks(cls, config):  # {{{
        """
        For each comparison grid, construct one task for computing the
        climatologies and one plotting task for each season.  The climatology
        task is a prerequisite of the plotting tasks, but the plotting tasks
        can run in parallel with one another.

        Parameters
        ----------
        config : MpasAnalysisConfigParser object
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        mocTask = cls(config=config)

        taskSuffix = 'MHT'
        seasons = ['ANN']

        variableList = ['timeMonthly_avg_meridionalHeatTransportLat',
                        'timeMonthly_avg_meridionalHeatTransportLatZ']

        climatologyTask = \
            MpasClimatology(config=config,
                            variableList=variableList,
                            taskSuffix=taskSuffix,
                            componentName='ocean',
                            seasons=seasons,
                            tags=['climatology'])

        # add climatologyTask as a prerequisite of the MOC task so
        # plotting won't happen until we have the required
        # climatologies
        if mocTask.prerequisiteTasks is None:
            mocTask.prerequisiteTasks = [climatologyTask.taskName]
        else:
            mocTask.prerequisiteTasks.append(climatologyTask.taskName)
        # We want to have access to some information from the
        # climatologyTask (namely, we need a way to find out what the
        # names of the climatology files are that it created), so we'll
        # keep a reference to it handy.
        mocTask.climatologyTask = climatologyTask

        tasks = [climatologyTask, mocTask]
        return tasks  # }}}

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
        super(MeridionalHeatTransport, self).__init__(
            config=config,
            taskName='meridionalHeatTransport',
            componentName='ocean',
            tags=['climatology'])

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
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(MeridionalHeatTransport, self).setup_and_check()

        config = self.config
        climatologyTask = self.climatologyTask

        self.check_analysis_enabled(
            analysisOptionName='config_am_meridionalheattransport_enable',
            raiseException=True)

        # call setup_and_check() on the climatology task because it will make
        # sure the start and end year are set and correct.  (In parallel mode,
        # this copy of the climatologyTask is different from the one that will
        # actually have run to completion before this task gets run.)
        climatologyTask.setup_and_check()

        self.startDate = climatologyTask.startDate
        self.endDate = climatologyTask.endDate
        self.startYear = climatologyTask.startYear
        self.endYear = climatologyTask.endYear

        # Later, we will read in depth and MHT latitude points
        # from mpaso.hist.am.meridionalHeatTransport.*.nc
        mhtFiles = self.historyStreams.readpath(
                'meridionalHeatTransportOutput')
        if len(mhtFiles) == 0:
            raise IOError('No MPAS-O MHT history file found: need at least '
                          'one ')

        self.mhtFile = mhtFiles[0]

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

        # }}}

    def run(self):  # {{{
        """
        Process MHT analysis member data if available.
        Plots MHT as:
           1D function of latitude
           2D function of latitude and depth

        Authors
        -------
        Mark Petersen, Milena Veneziani, Xylar Asay-Davis
        """
        print "\nPlotting meridional heat transport (MHT)..."

        config = self.config

        # Read in depth and MHT latitude points
        # Latitude is from binBoundaryMerHeatTrans written in
        #  mpaso.hist.am.meridionalHeatTransport.*.nc
        # Depth is from refZMid, also in
        # mpaso.hist.am.meridionalHeatTransport.*.nc

        print '  Read in depth and latitude...'
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

        print '\n  Compute and plot global meridional heat transport'

        print '   Load data...'

        # use the climatologyTask to get the right file name for the
        # computed climatology
        climatologyFileName = self.climatologyTask.get_ncclimo_file_name(
                season='ANN', stage='unmasked')

        annualClimatology = xr.open_dataset(climatologyFileName)

        # **** Plot MHT ****
        # Define plotting variables
        mainRunName = config.get('runs', 'mainRunName')
        xLimGlobal = config.getExpression(self.sectionName, 'xLimGlobal')
        depthLimGlobal = config.getExpression(self.sectionName,
                                              'depthLimGlobal')

        print '   Plot global MHT...'
        # Plot 1D MHT (zonally averaged, depth integrated)
        x = binBoundaryMerHeatTrans
        daY = annualClimatology.timeMonthly_avg_meridionalHeatTransportLat
        y = daY.values[0, :]
        xLabel = 'latitude [deg]'
        yLabel = 'meridional heat transport [PW]'
        title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
                 self.startYear, self.endYear, mainRunName)
        figureName = '{}/mht_{}_years{:04d}-{:04d}.png'.format(
                      self.plotsDirectory, mainRunName,
                      self.startYear, self.endYear)
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

        # Plot 2D MHT (zonally integrated)

        # normalize 2D MHT by layer thickness
        daMHTLatZ = \
            annualClimatology.timeMonthly_avg_meridionalHeatTransportLatZ
        MHTLatZ = daMHTLatZ.values[0, :, :].T
        for k in range(nVertLevels):
            MHTLatZ[k, :] = MHTLatZ[k, :]/refLayerThickness[k]

        x = binBoundaryMerHeatTrans
        y = refZMid
        z = MHTLatZ
        xLabel = 'latitude [deg]'
        yLabel = 'depth [m]'
        title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
                 self.startYear, self.endYear, mainRunName)
        figureName = '{}/mhtZ_{}_years{:04d}-{:04d}.png'.format(
                      self.plotsDirectory, mainRunName,
                      self.startYear, self.endYear)
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
                              invertYAxis=False)
        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
