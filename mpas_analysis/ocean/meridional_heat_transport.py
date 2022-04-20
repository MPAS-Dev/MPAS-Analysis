# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

import xarray as xr
import numpy as np
import os

from mpas_analysis.shared.plot import plot_vertical_section, plot_1D, savefig

from mpas_analysis.shared.io.utility import make_directories, build_obs_path
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.html import write_image_xml
from mpas_analysis.shared.climatology.climatology import \
    get_climatology_op_directory


class MeridionalHeatTransport(AnalysisTask):
    """
    Plot meridional heat transport from the analysis member output.

    Attributes
    ----------

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    controlconfig : mpas_tools.config.MpasConfigParser
        Configuration options for a control run (if any)
    """

    # Authors
    # -------
    # Mark Petersen, Milena Veneziani, Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(MeridionalHeatTransport, self).__init__(
            config=config,
            taskName='meridionalHeatTransport',
            componentName='ocean',
            tags=['climatology', 'publicObs'])

        self.mpasClimatologyTask = mpasClimatologyTask
        self.run_after(mpasClimatologyTask)

        self.controlConfig = controlConfig

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Mark Petersen, Milena Veneziani, Xylar Asay-Davis

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

        self.sectionName = 'meridionalHeatTransport'

        # Read in obs file information
        compareWithObs = config.getboolean(self.sectionName,
                                           'compareWithObservations')
        self.observationsFile = None
        if compareWithObs:
            observationsDirectory = build_obs_path(
                config, 'ocean', 'mhtSubdirectory')
            observationsFile = config.get(self.sectionName, 'observationData')
            observationsFile = '{}/{}'.format(observationsDirectory,
                                              observationsFile)
            if os.path.exists(observationsFile):
                self.observationsFile = observationsFile
            else:
                print('Warning: No MHT observations file found: skip plotting '
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

    def run_task(self):
        """
        Process MHT analysis member data if available.
        Plots MHT as:
           1D function of latitude
           2D function of latitude and depth
        """
        # Authors
        # -------
        # Mark Petersen, Milena Veneziani, Xylar Asay-Davis

        self.logger.info("\nPlotting meridional heat transport (MHT)...")

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')

        depthLimGlobal = config.getexpression(self.sectionName,
                                              'depthLimGlobal')
        xLimGlobal = config.getexpression(self.sectionName, 'xLimGlobal')
        movingAveragePoints = config.getint('meridionalHeatTransport',
                                            'movingAveragePoints')

        outputDirectory = get_climatology_op_directory(config)

        make_directories(outputDirectory)

        outFileName = \
            '{}/meridionalHeatTransport_years{:04d}-{:04d}.nc'.format(
                outputDirectory, self.startYear, self.endYear)

        if os.path.exists(outFileName):
            self.logger.info('  Reading results from previous analysis run...')
            annualClimatology = xr.open_dataset(outFileName)
            refZMid = annualClimatology.refZMid
            binBoundaryMerHeatTrans = \
                annualClimatology.binBoundaryMerHeatTrans
        else:

            # Read in depth and MHT latitude points
            # Latitude is from binBoundaryMerHeatTrans
            try:
                restartFileName = self.runStreams.readpath('restart')[0]
            except ValueError:
                raise IOError('No MPAS-O restart file found: need at least '
                              'one for MHT calcuation')

            with xr.open_dataset(restartFileName) as dsRestart:
                refBottomDepth = dsRestart.refBottomDepth

            nVertLevels = refBottomDepth.sizes['nVertLevels']
            refLayerThickness = np.zeros(nVertLevels)
            refLayerThickness[0] = refBottomDepth[0]
            refLayerThickness[1:nVertLevels] = \
                refBottomDepth[1:nVertLevels] - \
                refBottomDepth[0:nVertLevels - 1]

            refLayerThickness = xr.DataArray(dims='nVertLevels',
                                             data=refLayerThickness)

            refZMid = -refBottomDepth + 0.5 * refLayerThickness

            binBoundaryMerHeatTrans = None
            # first try timeSeriesStatsMonthly for bin boundaries, then try
            # meridionalHeatTransport stream as a backup option
            for streamName in ['timeSeriesStatsMonthlyOutput',
                               'meridionalHeatTransportOutput']:
                try:
                    inputFile = self.historyStreams.readpath(streamName)[0]
                except ValueError:
                    raise IOError('At least one file from stream {} is needed '
                                  'to compute MHT'.format(streamName))

                with xr.open_dataset(inputFile) as ds:
                    if 'binBoundaryMerHeatTrans' in ds.data_vars:
                        binBoundaryMerHeatTrans = \
                            ds.binBoundaryMerHeatTrans
                        break

            if binBoundaryMerHeatTrans is None:
                raise ValueError('Could not find binBoundaryMerHeatTrans in '
                                 'either timeSeriesStatsMonthlyOutput or '
                                 'meridionalHeatTransportOutput streams')

            binBoundaryMerHeatTrans = np.rad2deg(binBoundaryMerHeatTrans)

            ###################################################################
            # Mark P Note: Currently only supports global MHT.
            # Need to add variables merHeatTransLatRegion and
            # merHeatTransLatZRegion
            # These are not computed by default in ACME right now.
            # Then we will need to add another section for regions with a loop
            # over number of regions.
            ###################################################################

            self.logger.info('\n   Plotting global meridional heat transport')

            self.logger.info('   Load data...')

            climatologyFileName = self.mpasClimatologyTask.get_file_name(
                season='ANN')

            variableList = ['timeMonthly_avg_meridionalHeatTransportLat',
                            'timeMonthly_avg_meridionalHeatTransportLatZ']

            annualClimatology = xr.open_dataset(climatologyFileName)
            annualClimatology = annualClimatology[variableList]
            annualClimatology = annualClimatology.rename(
                {'timeMonthly_avg_meridionalHeatTransportLat':
                     'meridionalHeatTransportLat',
                 'timeMonthly_avg_meridionalHeatTransportLatZ':
                     'meridionalHeatTransportLatZ'})
            if 'Time' in annualClimatology.dims:
                annualClimatology = annualClimatology.isel(Time=0)

            annualClimatology.coords['refZMid'] = refZMid
            annualClimatology.coords['binBoundaryMerHeatTrans'] = \
                binBoundaryMerHeatTrans

            if config.getboolean(self.sectionName, 'plotVerticalSection'):
                # normalize 2D MHT by layer thickness
                annualClimatology['meridionalHeatTransportLatZ'] /= \
                    refLayerThickness

            write_netcdf(annualClimatology, outFileName)

        # **** Plot MHT ****
        maxTitleLength = 70
        self.logger.info('   Plot global MHT...')
        # Plot 1D MHT (zonally averaged, depth integrated)
        x = binBoundaryMerHeatTrans
        y = annualClimatology.meridionalHeatTransportLat
        xLabel = 'latitude (deg)'
        yLabel = 'meridional heat transport (PW)'

        title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
            self.startYear, self.endYear, mainRunName)
        filePrefix = self.filePrefixes['mht']
        figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        lineColors = ['k']
        lineWidths = [1.6]
        legendText = [mainRunName]
        xArrays = [x]
        fieldArrays = [y]
        errArrays = [None]
        if self.observationsFile is not None:
            # Load in observations
            dsObs = xr.open_dataset(self.observationsFile)
            xObs = dsObs.LATITUDE
            ncepGlobal = dsObs.GLOBALNCEP_ADJUSTED
            ncepErrGlobal = dsObs.GLOBALNCEP_ERR
            ecmwfGlobal = dsObs.GLOBALECMWF_ADJUSTED
            ecmwfErrGlobal = dsObs.GLOBALECMWF_ERR

            lineColors.extend(['b', 'g'])
            lineWidths.extend([1.2, 1.2])
            legendText.extend(['Trenberth and Caron - NCEP',
                               'Trenberth and Caron - ECMWF'])
            xArrays.extend([xObs, xObs])
            fieldArrays.extend([ncepGlobal, ecmwfGlobal])
            errArrays.extend([ncepErrGlobal, ecmwfErrGlobal])

        if self.controlConfig is not None:
            controlStartYear = self.controlConfig.getint('climatology',
                                                         'startYear')
            controlEndYear = self.controlConfig.getint('climatology',
                                                       'endYear')
            controlDirectory = get_climatology_op_directory(self.controlConfig)

            controlFileName = \
                '{}/meridionalHeatTransport_years{:04d}-{:04d}.nc'.format(
                    controlDirectory, controlStartYear, controlEndYear)

            dsControl = xr.open_dataset(controlFileName)
            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            lineColors.append('r')
            lineWidths.append(1.2)
            legendText.append(controlRunName)
            xArrays.append(dsControl.binBoundaryMerHeatTrans)
            fieldArrays.append(dsControl.meridionalHeatTransportLat)
            errArrays.append(None)

        if len(legendText) == 1:
            # no need for a legend
            legendText = [None]

        plot_1D(config, xArrays, fieldArrays, errArrays,
                lineColors=lineColors, lineWidths=lineWidths,
                legendText=legendText, title=title, xlabel=xLabel,
                ylabel=yLabel, fileout=figureName, xLim=xLimGlobal,
                maxTitleLength=maxTitleLength)

        self._write_xml(filePrefix)

        if config.getboolean(self.sectionName, 'plotVerticalSection'):
            # Plot 2D MHT (zonally integrated)

            x = binBoundaryMerHeatTrans
            y = refZMid
            z = annualClimatology.meridionalHeatTransportLatZ
            xLabel = 'latitude (deg)'
            yLabel = 'depth (m)'
            title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
                self.startYear, self.endYear, mainRunName)
            filePrefix = self.filePrefixes['mhtZ']
            outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
            colorbarLabel = '(PW/m)'
            plot_vertical_section(config, z, self.sectionName, xCoords=x,
                                  zCoord=y, suffix='',
                                  colorbarLabel=colorbarLabel,
                                  title=title, xlabels=xLabel, ylabel=yLabel,
                                  xLim=xLimGlobal,
                                  yLim=depthLimGlobal, invertYAxis=False,
                                  movingAveragePoints=movingAveragePoints,
                                  maxTitleLength=maxTitleLength)

            savefig(outFileName, config)

            self._write_xml(filePrefix)

    def _write_xml(self, filePrefix):
        caption = 'Meridional Heat Transport'
        write_image_xml(
            config=self.config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Meridional Heat Transport',
            groupLink='mht',
            imageDescription=caption,
            imageCaption=caption)
