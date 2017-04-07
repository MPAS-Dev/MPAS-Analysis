"""
Plot meridional heat transport from the analysis member output.

Authors
-------
Mark Petersen, Milena Veneziani

Last Modified
-------------
04/21/2017
"""

import xarray as xr
import numpy as np
import netCDF4
import os
import warnings

from ..shared.constants.constants import rad_to_deg, monthDictionary
from ..shared.plot.plotting import plot_vertical_section,\
    setup_colormap, plot_1D

from ..shared.io.utility import build_config_full_path, make_directories

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.analysis_task import setup_task, check_analysis_enabled
from ..shared.climatology.climatology import cache_climatologies


def meridional_heat_transport(config):  # {{{
    """
    Process MHT analysis member data if available.
    Plots MHT as:
       1D function of latitude
       2D function of latitude and depth

    Parameters
    ----------
    config : instance of MpasAnalysisConfigParser
        configuration options used to customize the analysis task

    Authors
    -------
    Mark Petersen, Milena Veneziani

    Last Modified
    -------------
    04/21/2017
    """

    # **** Initial settings ****

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    check_analysis_enabled(
        namelist=namelist,
        analysisOptionName='config_am_timeseriesstatsmonthly_enable',
        namelistMap=namelistMap,
        raiseException=True)
    mhtAnalysisMemberEnabled = check_analysis_enabled(
        namelist=namelist,
        analysisOptionName='config_am_meridionalheattransport_enable',
        namelistMap=namelistMap,
        raiseException=False)
    if not mhtAnalysisMemberEnabled:
        message = '*** MHT Analysis Member is off ***\n' \
                  '***     Unable to plot MHT     ***\n'
        warnings.warn(message)
        return

    # Get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    #   First a list necessary for theMHT climatology
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    startDateClimo = config.get('climatology', 'startDate')
    endDateClimo = config.get('climatology', 'endDate')
    inputFilesClimo = historyStreams.readpath(streamName,
                                              startDate=startDateClimo,
                                              endDate=endDateClimo,
                                              calendar=calendar)
    simulationStartTime = get_simulation_start_time(runStreams)

    print '\n  List of files for climatologies:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(inputFilesClimo[0]),
              os.path.basename(inputFilesClimo[-1]))

    startYearClimo = config.getint('climatology', 'startYear')
    endYearClimo = config.getint('climatology', 'endYear')

    sectionName = 'meridionalHeatTransport'

    # Read in obs file information
    compareWithObs = config.getboolean(sectionName, 'compareWithObservations')
    if compareWithObs is True:
        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'mhtSubdirectory')
        observationsFile = config.get(sectionName, 'observationData')
        observationsFile = '{}/{}'.format(observationsDirectory,
                                          observationsFile)
        if not os.path.exists(observationsFile):
            warnings('No MHT observations file found: skip plotting obs')
            compareWithObs = False

    # Read in depth and MHT latitude points
    # Latitude is from binBoundaryMerHeatTrans written in
    #  mpaso.hist.am.meridionalHeatTransport.*.nc
    # Depth is from refZMid, also in mpaso.hist.am.meridionalHeatTransport.*.nc
    try:
        mhtFile = historyStreams.readpath('meridionalHeatTransportOutput')[0]
    except ValueError:
        raise IOError('No MPAS-O MHT history file found: need at least one ')

    print '  Read in depth and latitude...'
    ncFile = netCDF4.Dataset(mhtFile, mode='r')
    # reference depth [m]
    refZMid = ncFile.variables['refZMid'][:]
    refBottomDepth = ncFile.variables['refBottomDepth'][:]
    binBoundaryMerHeatTrans = ncFile.variables['binBoundaryMerHeatTrans'][:]
    binBoundaryMerHeatTrans = binBoundaryMerHeatTrans*rad_to_deg
    ncFile.close()

    nVertLevels = len(refBottomDepth)
    refLayerThickness = np.zeros(nVertLevels)
    refLayerThickness[0] = refBottomDepth[0]
    refLayerThickness[1:nVertLevels] = (refBottomDepth[1:nVertLevels] -
                                        refBottomDepth[0:nVertLevels-1])

    ########################################################################
    # Mark P Note: Currently only supports global MHT.
    # Need to add variables merHeatTransLatRegion and merHeatTransLatZRegion
    # these are not computed by default in ACME right now.
    # Then we will need to add another section for regions with a loop
    # over number of regions.
    ########################################################################
    variableList = ['avgMeridionalHeatTransportLat',
                    'avgMeridionalHeatTransportLatZ']

    print '\n  Compute and plot global meridional heat transport'

    outputDirectory = build_config_full_path(config, 'output',
                                             'mpasClimatologySubdirectory')

    make_directories(outputDirectory)

    print '   Load data...'
    ds = open_multifile_dataset(
        fileNames=inputFilesClimo,
        calendar=calendar,
        config=config,
        simulationStartTime=simulationStartTime,
        timeVariableName='Time',
        variableList=variableList,
        variableMap=variableMap,
        startDate=startDateClimo,
        endDate=endDateClimo)

    # Compute annual climatology
    cachePrefix = '{}/meridionalHeatTransport'.format(outputDirectory)
    annualClimatology = cache_climatologies(ds, monthDictionary['ANN'],
                                            config, cachePrefix,
                                            calendar, printProgress=True)

    # **** Plot MHT ****
    # Define plotting variables
    mainRunName = config.get('runs', 'mainRunName')
    xLimGlobal = config.getExpression(sectionName, 'xLimGlobal')
    depthLimGlobal = config.getExpression(sectionName, 'depthLimGlobal')

    print '   Plot global MHT...'
    # Plot 1D MHT (zonally averaged, depth integrated)
    x = binBoundaryMerHeatTrans
    y = annualClimatology.avgMeridionalHeatTransportLat
    xLabel = 'latitude [deg]'
    yLabel = 'meridional heat transport [PW]'
    title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
             startYearClimo, endYearClimo, mainRunName)
    figureName = '{}/mht_{}_years{:04d}-{:04d}.png'.format(
                  plotsDirectory, mainRunName,
                  startYearClimo, endYearClimo)
    if compareWithObs is True:
        # Load in observations
        dsObs = xr.open_dataset(observationsFile)
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
    MHTLatZ = annualClimatology.avgMeridionalHeatTransportLatZ.values.T[:,:]
    for k in range(nVertLevels):
    	MHTLatZ[k,:] = MHTLatZ[k,:]/refLayerThickness[k]

    x = binBoundaryMerHeatTrans
    y = refZMid
    z = MHTLatZ
    xLabel = 'latitude [deg]'
    yLabel = 'depth [m]'
    title = 'Global MHT (ANN, years {:04d}-{:04d})\n {}'.format(
             startYearClimo, endYearClimo, mainRunName)
    figureName = '{}/mhtZ_{}_years{:04d}-{:04d}.png'.format(
                  plotsDirectory, mainRunName,
                  startYearClimo, endYearClimo)
    colorbarLabel = '[PW/m]'
    contourLevels = config.getExpression(sectionName,
                                         'contourLevelsGlobal',
                                         usenumpyfunc=True)
    (colormapName, colorbarLevels) = setup_colormap(config, sectionName,
                                                    suffix='Global')
    plot_vertical_section(config, x, y, z,
                          colormapName, colorbarLevels,
                          contourLevels, colorbarLabel,
                          title, xLabel, yLabel, figureName,
                          xLim=xLimGlobal, yLim=depthLimGlobal,
                          invertYAxis=False)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
