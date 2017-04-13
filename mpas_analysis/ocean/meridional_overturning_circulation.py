"""
Computation and plotting of model meridional overturning circulation.
Will eventually support:
  * MOC streamfunction, post-processed (currently supported)
  * MOC streamfunction, from MOC analysis member
  * MOC time series (max value at 24.5N), post-processed
  * MOC time series (max value at 24.5N), from MOC analysis member

Authors
-------
Milena Veneziani, Mark Petersen, Phillip Wolfram, Xylar Asay-Davis

Last Modified
-------------
04/08/2017
"""

import xarray as xr
import numpy as np
import netCDF4
import os
from functools import partial

from ..shared.constants.constants import m3ps_to_Sv, rad_to_deg, \
    monthDictionary
from ..shared.plot.plotting import plot_vertical_section,\
    timeseries_analysis_plot, setup_colormap

from ..shared.io.utility import build_config_full_path, make_directories

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    days_to_datetime

from ..shared.analysis_task import setup_task, check_analysis_enabled
from ..shared.climatology.climatology import update_start_end_year, \
    cache_climatologies
from ..shared.time_series import time_series


def moc_streamfunction(config):  # {{{
    """
    Process MOC analysis member data if available, or compute MOC at
    post-processing if not. Plots streamfunction climatolgoical sections
    as well as time series of max Atlantic MOC at 26.5N (latitude of
    RAPID MOC Array).

    Parameters
    ----------
    config : instance of MpasAnalysisConfigParser
        configuration options used to customize the analysis task

    streamMap : dict, optional
        a dictionary of MPAS-O stream names that map to
        their mpas_analysis counterparts.

    variableMap : dict, optional
        a dictionary of MPAS-O variable names that map
        to their mpas_analysis counterparts.

    Authors
    -------
    Milena Veneziani, Mark Petersen, Phillip J. Wolfram, Xylar Asay-Davis

    Last Modified
    -------------
    04/08/2017
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

    # Get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    #   First a list necessary for the streamfunctionMOC climatology
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
    #   Create dictionary to store Climo related variables
    dictClimo = {'inputFilesClimo': inputFilesClimo,
                 'startDateClimo': startDateClimo,
                 'endDateClimo':  endDateClimo,
                 'startYearClimo': startYearClimo,
                 'endYearClimo': endYearClimo}

    #   Then a list necessary for the streamfunctionMOC Atlantic timeseries
    startDateTseries = config.get('timeSeries', 'startDate')
    endDateTseries = config.get('timeSeries', 'endDate')
    inputFilesTseries = historyStreams.readpath(streamName,
                                                startDate=startDateTseries,
                                                endDate=endDateTseries,
                                                calendar=calendar)

    print '\n  List of files for time series:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(inputFilesTseries[0]),
              os.path.basename(inputFilesTseries[-1]))

    startYearTseries = config.getint('timeSeries', 'startYear')
    endYearTseries = config.getint('timeSeries', 'endYear')
    #   Create dictionary to store Tseries related variables
    dictTseries = {'inputFilesTseries': inputFilesTseries,
                   'startDateTseries': startDateTseries,
                   'endDateTseries':  endDateTseries,
                   'startYearTseries': startYearTseries,
                   'endYearTseries': endYearTseries}

    sectionName = 'streamfunctionMOC'
    regionNames = config.getExpression(sectionName, 'regionNames')

    # **** Compute MOC ****
    mocAnalysisMemberEnabled = check_analysis_enabled(
        namelist=namelist,
        analysisOptionName='config_am_mocstreamfunction_enable',
        namelistMap=namelistMap,
        raiseException=False)

    # Check whether MOC Analysis Member is enabled
    if mocAnalysisMemberEnabled:
        # Add a moc_analisysMember_processing
        print '*** MOC Analysis Member is on ***'
        # (mocDictClimo, mocDictTseries) = _compute_moc_analysismember(config,
        #     streams, calendar, sectionName, dictClimo, dictTseries)
    else:
        _cache_velocity_climatologies(config, sectionName,
                                      startDateClimo, endDateClimo,
                                      inputFilesClimo, simulationStartTime,
                                      variableMap, calendar)

        # update the start and end year in case they have changed
        dictClimo['startYearClimo'] = config.getint('climatology', 'startYear')
        dictClimo['endYearClimo'] = config.getint('climatology', 'endYear')

        mocDictClimo, dictRegion = _compute_moc_climo_postprocess(
            config, runStreams, variableMap, calendar, sectionName,
            regionNames, dictClimo)
        dsMOCTimeSeries = _compute_moc_time_series_postprocess(
            config, runStreams, variableMap, calendar, sectionName,
            regionNames, dictTseries, mocDictClimo, dictRegion)

    # **** Plot MOC ****
    # Define plotting variables
    mainRunName = config.get('runs', 'mainRunName')
    movingAveragePoints = config.getint(sectionName, 'movingAveragePoints')
    colorbarLabel = '[Sv]'
    xLabel = 'latitude [deg]'
    yLabel = 'depth [m]'

    for region in regionNames:
        print '   Plot climatological {} MOC...'.format(region)
        title = '{} MOC (ANN, years {:04d}-{:04d})\n {}'.format(
                 region, dictClimo['startYearClimo'],
                 dictClimo['endYearClimo'],
                 mainRunName)
        figureName = '{}/moc{}_{}_years{:04d}-{:04d}.png'.format(
                      plotsDirectory, region, mainRunName,
                      dictClimo['startYearClimo'], dictClimo['endYearClimo'])
        contourLevels = config.getExpression(sectionName,
                                             'contourLevels{}'.format(region),
                                             usenumpyfunc=True)
        (colormapName, colorbarLevels) = setup_colormap(config, sectionName,
                                                        suffix=region)

        x = mocDictClimo['lat{}'.format(region)]['data']
        y = mocDictClimo['depth']['data']
        z = mocDictClimo['moc{}'.format(region)]['data']
        plot_vertical_section(config, x, y, z, colormapName, colorbarLevels,
                              contourLevels, colorbarLabel, title,
                              xLabel, yLabel, figureName)

    # Plot time series
    print '   Plot time series of max Atlantic MOC at 26.5N...'
    xLabel = 'Time [years]'
    yLabel = '[Sv]'
    title = 'Max Atlantic MOC at $26.5^\circ$N\n {}'.format(mainRunName)
    figureName = '{}/mocTimeseries_{}.png'.format(plotsDirectory,
                                                  mainRunName)

    timeseries_analysis_plot(config, [dsMOCTimeSeries.mocAtlantic26],
                             movingAveragePoints, title,
                             xLabel, yLabel, figureName,
                             lineStyles=['k-'], lineWidths=[1.5],
                             calendar=calendar)
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
    latCell = ncFile.variables['latCell'][:]
    latCell = latCell * rad_to_deg  # convert to degree
    ncFile.close()
    nVertLevels = len(refBottomDepth)
    refTopDepth = np.zeros(nVertLevels+1)
    refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]
    refLayerThickness = np.zeros(nVertLevels)
    refLayerThickness[0] = refBottomDepth[0]
    refLayerThickness[1:nVertLevels] = (refBottomDepth[1:nVertLevels] -
                                        refBottomDepth[0:nVertLevels-1])

    return dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        refTopDepth, refLayerThickness  # }}}


def _cache_velocity_climatologies(config, sectionName,
                                  startDateClimo, endDateClimo,
                                  inputFilesClimo, simulationStartTime,
                                  variableMap, calendar):  # {{{
    '''compute yearly velocity climatologies and cache them'''

    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']

    outputDirectory = build_config_full_path(config, 'output',
                                             'mpasClimatologySubdirectory')

    make_directories(outputDirectory)

    chunking = config.getExpression(sectionName, 'maxChunkSize')
    ds = open_multifile_dataset(fileNames=inputFilesClimo,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=startDateClimo,
                                endDate=endDateClimo,
                                chunking=chunking)

    # update the start and end year in config based on the real extend of ds
    update_start_end_year(ds, config, calendar)

    cachePrefix = '{}/meanVelocity'.format(outputDirectory)

    # compute and cache the velocity climatology
    cache_climatologies(ds, monthDictionary['ANN'],
                        config, cachePrefix, calendar,
                        printProgress=True)
    # }}}


def _compute_moc_climo_postprocess(config, runStreams, variableMap, calendar,
                                   sectionName, regionNames, dictClimo):  # {{{

    '''compute mean MOC streamfunction as a post-process'''

    dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        refTopDepth, refLayerThickness = _load_mesh(runStreams)

    # Load basin region related variables and save them to dictionary
    # NB: The following will need to change with new regional mapping files
    regionMaskFiles = config.get(sectionName, 'regionMaskFiles')
    if not os.path.exists(regionMaskFiles):
        raise IOError('Regional masking file for MOC calculation '
                      'does not exist')
    iRegion = 0
    for region in regionNames:
        print '\n  Reading region and transect mask for {}...'.format(region)
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
        dictRegion = {
            'ind{}'.format(region): indRegion,
            '{}CellMask'.format(region): regionCellMask,
            'maxEdgesInTransect{}'.format(region): maxEdgesInTransect,
            'transectEdgeMaskSigns{}'.format(region): transectEdgeMaskSigns,
            'transectEdgeGlobalIDs{}'.format(region): transectEdgeGlobalIDs}
    # Add Global regionCellMask=1 everywhere to make the algorithm
    # for the global moc similar to that of the regional moc
    dictRegion['GlobalCellMask'] = np.ones(np.size(latCell))
    regionNames[:0] = ['Global']

    # Compute and plot annual climatology of MOC streamfunction
    print '\n  Compute and/or plot post-processed MOC climatological '\
          'streamfunction...'
    outputDirectory = build_config_full_path(config, 'output',
                                             'mpasClimatologySubdirectory')

    make_directories(outputDirectory)

    outputFileClimo = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
                       outputDirectory, dictClimo['startYearClimo'],
                       dictClimo['endYearClimo'])
    if not os.path.exists(outputFileClimo):
        print '   Load data...'

        velClimoFile = '{}/meanVelocity_years{:04d}-{:04d}.nc'.format(
                       outputDirectory, dictClimo['startYearClimo'],
                       dictClimo['endYearClimo'])

        annualClimatology = xr.open_dataset(velClimoFile)

        # Convert to numpy arrays
        # (can result in a memory error for large array size)
        horizontalVel = annualClimatology.avgNormalVelocity.values
        verticalVel = annualClimatology.avgVertVelocityTop.values
        velArea = verticalVel * areaCell[:, np.newaxis]

        # Create dictionary for MOC climatology (NB: need this form
        # in order to convert it to xarray dataset later in the script)
        mocDictClimo = {'depth': {'dims': ('nz'), 'data': refTopDepth}}
        for region in regionNames:
            print '   Compute {} MOC...'.format(region)
            print '    Compute transport through region southern transect...'
            if region == 'Global':
                transportZ = np.zeros(nVertLevels)
            else:
                maxEdgesInTransect = \
                    dictRegion['maxEdgesInTransect{}'.format(region)]
                transectEdgeGlobalIDs = \
                    dictRegion['transectEdgeGlobalIDs{}'.format(region)]
                transectEdgeMaskSigns = \
                    dictRegion['transectEdgeMaskSigns{}'.format(region)]
                transportZ = _compute_transport(maxEdgesInTransect,
                                                transectEdgeGlobalIDs,
                                                transectEdgeMaskSigns,
                                                nVertLevels, dvEdge,
                                                refLayerThickness,
                                                horizontalVel)

            regionCellMask = dictRegion['{}CellMask'.format(region)]
            latBinSize = config.getExpression(sectionName,
                                              'latBinSize{}'.format(region))
            if region == 'Global':
                latBins = np.arange(-90.0, 90.1, latBinSize)
            else:
                indRegion = dictRegion['ind{}'.format(region)]
                latBins = latCell[indRegion]
                latBins = np.arange(np.amin(latBins),
                                    np.amax(latBins)+latBinSize,
                                    latBinSize)
            mocTop = _compute_moc(latBins, nVertLevels, latCell,
                                  regionCellMask, transportZ, velArea)

            # Store computed MOC to dictionary
            mocDictClimo['lat{}'.format(region)] = {
                'dims': ('nx{}'.format(region)), 'data': latBins}
            mocDictClimo['moc{}'.format(region)] = {
                'dims': ('nz', 'nx{}'.format(region)), 'data': mocTop}

        # Save to file
        print '   Save global and regional MOC to file...'
        ncFile = netCDF4.Dataset(outputFileClimo, mode='w')
        # create dimensions
        ncFile.createDimension('nz', len(refTopDepth))
        for region in regionNames:
            latBins = mocDictClimo['lat{}'.format(region)]['data']
            mocTop = mocDictClimo['moc{}'.format(region)]['data']
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
        refTopDepth = ncFile.variables['depth'][:]
        mocDictClimo = {'depth': {'dims': ('nz'), 'data': refTopDepth}}
        for region in regionNames:
            latBins = ncFile.variables['lat{}'.format(region)][:]
            mocTop = ncFile.variables['moc{}'.format(region)][:, :]
            mocDictClimo['lat{}'.format(region)] = {
                'dims': ('nx{}'.format(region)), 'data': latBins}
            mocDictClimo['moc{}'.format(region)] = {
                'dims': ('nz', 'nx{}'.format(region)), 'data': mocTop}
        ncFile.close()
    return mocDictClimo, dictRegion  # }}}


def _compute_moc_time_series_postprocess(config, runStreams, variableMap,
                                         calendar, sectionName, regionNames,
                                         dictTseries, mocDictClimo,
                                         dictRegion):  # {{{
    '''compute MOC time series as a post-process'''

    # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
    print '\n  Compute and/or plot post-processed Atlantic MOC '\
          'time series...'
    print '   Load data...'

    simulationStartTime = get_simulation_start_time(runStreams)
    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']

    dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        refTopDepth, refLayerThickness = _load_mesh(runStreams)

    chunking = config.getExpression(sectionName, 'maxChunkSize')
    ds = open_multifile_dataset(fileNames=dictTseries['inputFilesTseries'],
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=dictTseries['startDateTseries'],
                                endDate=dictTseries['endDateTseries'],
                                chunking=chunking)
    latAtlantic = mocDictClimo['latAtlantic']['data']
    dLat = latAtlantic - 26.5
    indlat26 = np.where(dLat == np.amin(np.abs(dLat)))

    maxEdgesInTransect = dictRegion['maxEdgesInTransectAtlantic']
    transectEdgeGlobalIDs = dictRegion['transectEdgeGlobalIDsAtlantic']
    transectEdgeMaskSigns = dictRegion['transectEdgeMaskSignsAtlantic']
    regionCellMask = dictRegion['AtlanticCellMask']

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
    comp_moc_part = partial(_compute_moc_time_series_part, ds,
                            calendar, areaCell, latCell, indlat26,
                            maxEdgesInTransect, transectEdgeGlobalIDs,
                            transectEdgeMaskSigns,  nVertLevels, dvEdge,
                            refLayerThickness, latAtlantic, regionCellMask)

    dsMOCTimeSeries = time_series.cache_time_series(
        ds.Time.values,  comp_moc_part, outputFileTseries,
        calendar, yearsPerCacheUpdate=1,  printProgress=False)

    return dsMOCTimeSeries  # }}}


def _compute_moc_time_series_part(ds, calendar, areaCell, latCell, indlat26,
                                  maxEdgesInTransect, transectEdgeGlobalIDs,
                                  transectEdgeMaskSigns, nVertLevels, dvEdge,
                                  refLayerThickness, latAtlantic,
                                  regionCellMask, timeIndices, firstCall):
    # computes a subset of the MOC time series

    if firstCall:
        print '   Process and save time series'

    times = ds.Time[timeIndices].values
    mocRegion = np.zeros(timeIndices.shape)

    for localIndex, timeIndex in enumerate(timeIndices):
        time = times[localIndex]
        dsLocal = ds.isel(Time=timeIndex)
        date = days_to_datetime(time, calendar=calendar)

        print '     date: {:04d}-{:02d}'.format(date.year, date.month)

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


# def _compute_moc_analysismember(config):
#
#     return (mocDictClimo, mocDictTseries)


def _compute_transport(maxEdgesInTransect, transectEdgeGlobalIDs,
                       transectEdgeMaskSigns, nz, dvEdge, refLayerThickness,
                       horizontalVel):  # {{{

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

    mocTop = np.zeros([np.size(latBins), nz+1])
    mocTop[0, range(1, nz+1)] = transportZ.cumsum()
    for iLat in range(1, np.size(latBins)):
        indlat = np.logical_and(np.logical_and(
                     regionCellMask == 1, latCell >= latBins[iLat-1]),
                     latCell < latBins[iLat])
        mocTop[iLat, :] = mocTop[iLat-1, :] + velArea[indlat, :].sum(axis=0)
    # convert m^3/s to Sverdrup
    mocTop = mocTop * m3ps_to_Sv
    mocTop = mocTop.T
    return mocTop  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
