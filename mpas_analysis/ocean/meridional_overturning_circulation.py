"""
Computation and plotting of model meridional overturning circulation.
Will eventually support:
  * MOC streamfunction, post-processed (currently supported)
  * MOC streamfunction, from MOC analysis member
  * MOC time series (max value at 24.5N), post-processed
  * MOC time series (max value at 24.5N), from MOC analysis member

Authors
-------
Milena Veneziani, Mark Petersen, Phillip Wolfram

Last Modified
-------------
03/14/2017
"""

import xarray as xr
import numpy as np
import netCDF4
import os

from ..shared.constants.constants import m3ps_to_Sv, rad_to_deg
from ..shared.plot.plotting import plot_vertical_section,\
    timeseries_analysis_plot, setup_colormap

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    days_to_datetime

from ..shared.climatology import climatology


def moc_streamfunction(config, streamMap=None, variableMap=None):  # {{{
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
    Milena Veneziani, Mark Petersen, Phillip J. Wolfram

    Last Modified
    -------------
    03/14/2017
    """

    # **** Initial settings ****
    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    timeSeriesStatsMonthlyAnalysisMemberFlag = namelist.get(
        'config_am_timeseriesstatsmonthly_enable')
    if timeSeriesStatsMonthlyAnalysisMemberFlag == '.false.':
        raise RuntimeError('*** MPAS-Analysis relies on the existence of monthly\n'
                           '*** averaged files: make sure to enable the\n'
                           '*** timeSeriesStatsMonthly analysis member!')

    # Get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    #   First a list necessary for the streamfunctionMOC climatology
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    startDateClimo = config.get('climatology', 'startDate')
    endDateClimo = config.get('climatology', 'endDate')
    calendar = namelist.get('config_calendar_type')
    inputFilesClimo = streams.readpath(streamName,
                                       startDate=startDateClimo,
                                       endDate=endDateClimo,
                                       calendar=calendar)
    print '\n  List of files for climatologies:\n{} through {}'.format(
           inputFilesClimo[0], inputFilesClimo[-1])
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
    inputFilesTseries = streams.readpath(streamName,
                                         startDate=startDateTseries,
                                         endDate=endDateTseries,
                                         calendar=calendar)
    print '\n  List of files for timeSeries:\n{} through {}'.format(
           inputFilesTseries[0], inputFilesTseries[-1])
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
    mocAnalysisMemberFlag = namelist.get('config_am_mocstreamfunction_enable')
    # Check whether MOC Analysis Member is enabled
    if mocAnalysisMemberFlag == '.true.':
        # Add a moc_analisysMember_processing
        print '*** MOC Analysis Member is on ***'
        # (mocDictClimo, mocDictTseries) = _compute_moc_analysismember(config,
        #     streams, calendar, sectionName, dictClimo, dictTseries)
    else:
        (mocDictClimo, mocDictTseries) = _compute_moc_postprocess(
            config, streams, variableMap, calendar, sectionName, regionNames,
            dictClimo, dictTseries)

    # **** Plot MOC ****
    # Define plotting variables
    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')
    movingAveragePoints = config.getint(sectionName, 'movingAveragePoints')
    colorbarLabel = '[Sv]'
    xLabel = 'latitude [deg]'
    yLabel = 'depth [m]'

    for region in regionNames:
        print '   Plot climatological {} MOC...'.format(region)
        title = '{} MOC (ANN, years {:04d}-{:04d})\n {}'.format(
                 region, startYearClimo, endYearClimo, mainRunName)
        figureName = '{}/moc{}_{}_years{:04d}-{:04d}.png'.format(
                      plotsDirectory, region, mainRunName,
                      startYearClimo, endYearClimo)
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

    dsmoc = xr.Dataset.from_dict(mocDictTseries)
    timeseries_analysis_plot(config, [dsmoc.mocRapidArrayMax],
                             movingAveragePoints, title,
                             xLabel, yLabel, figureName,
                             lineStyles=['k-'], lineWidths=[1.5],
                             calendar=calendar)

    return  # }}}


def _compute_moc_postprocess(config, streams, variableMap, calendar,
                             sectionName, regionNames, dictClimo,
                             dictTseries):

    '''compute MOC streamfunction at postprocessing'''

    # Load mesh related variables
    try:
        restartFile = streams.readpath('restart')[0]
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

    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']

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
    simulationStartTime = get_simulation_start_time(streams)
    outputDirectory = buildConfigFullPath(config, 'output',
                                          'mpasClimatologySubdirectory')
    try:
        os.makedirs(outputDirectory)
    except OSError:
        pass
    outputFileClimo = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
                       outputDirectory, dictClimo['startYearClimo'],
                       dictClimo['endYearClimo'])
    if not os.path.exists(outputFileClimo):
        print '   Load data...'
        ds = open_multifile_dataset(
            fileNames=dictClimo['inputFilesClimo'],
            calendar=calendar,
            simulationStartTime=simulationStartTime,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=variableMap,
            startDate=dictClimo['startDateClimo'],
            endDate=dictClimo['endDateClimo'])

        # Compute annual climatology
        annualClimatology = climatology.compute_annual_climatology(
                            ds, calendar)

        horizontalVel = annualClimatology.avgNormalVelocity
        verticalVel = annualClimatology.avgVertVelocityTop
        # Convert to numpy
        horizontalVel = horizontalVel.values
        verticalVel = verticalVel.values
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

    # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
    print '\n  Compute and/or plot post-processed Atlantic MOC '\
          'time series...'
    print '   Load data...'
    ds = open_multifile_dataset(fileNames=dictTseries['inputFilesTseries'],
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=dictTseries['startDateTseries'],
                                endDate=dictTseries['endDateTseries'])
    latAtlantic = mocDictClimo['latAtlantic']['data']
    dLat = latAtlantic - 26.5
    indlat26 = np.where(dLat == np.amin(np.abs(dLat)))

    startYear = days_to_datetime(ds.Time.min()).year
    endDate = days_to_datetime(ds.Time.max())
    endYear = endDate.year
    if endDate.month < 12:
        endYear -= 1
    nyears = (endYear - startYear + 1)

    Time = np.zeros(nyears*12)
    mocAtlantic26 = np.zeros(nyears*12)
    maxEdgesInTransect = dictRegion['maxEdgesInTransectAtlantic']
    transectEdgeGlobalIDs = dictRegion['transectEdgeGlobalIDsAtlantic']
    transectEdgeMaskSigns = dictRegion['transectEdgeMaskSignsAtlantic']
    regionCellMask = dictRegion['AtlanticCellMask']

    outputDirectory = buildConfigFullPath(config, 'output',
                                          'timeseriesSubdirectory')
    try:
        os.makedirs(outputDirectory)
    except OSError:
        pass

    ncount = 0
    for ny in range(nyears):
        year = dictTseries['startYearTseries'] + ny
        print '   year = ', year
        outputFileTseries = '{}/mocTimeSeries_year{:04d}.nc'.format(
                             outputDirectory, year)
        if not os.path.exists(outputFileTseries):
            for nm in range(12):
                print '     month = ', nm + 1
                horizontalVel = ds.avgNormalVelocity[ncount, :, :].values
                verticalVel = ds.avgVertVelocityTop[ncount, :, :].values
                velArea = verticalVel * areaCell[:, np.newaxis]
                transportZ = _compute_transport(maxEdgesInTransect,
                                                transectEdgeGlobalIDs,
                                                transectEdgeMaskSigns,
                                                nVertLevels, dvEdge,
                                                refLayerThickness,
                                                horizontalVel)
                mocTop = _compute_moc(latAtlantic, nVertLevels, latCell,
                                      regionCellMask, transportZ, velArea)
                mocAtlantic26[ncount] = np.amax(mocTop[:, indlat26])
                ncount += 1

            Time[ncount-12:ncount] = ds.Time[ncount-12:ncount].values

            # Save to file
            print '   save to yearly file...'
            ncFile = netCDF4.Dataset(outputFileTseries, mode='w')
            # create dimension
            ncFile.createDimension('Time', 0)
            # create variables
            Timevar = ncFile.createVariable('Time', 'f4', ('Time',))
            Timevar.description = 'days since start of simulation'
            mocAtlantic26var = ncFile.createVariable('mocRapidArrayMax',
                                                     'f4', ('Time',))
            mocAtlantic26var.description = 'Max MOC Atlantic streamfunction'\
                                           ' nearest to RAPID Array latitude'\
                                           ' (26.5N)'
            mocAtlantic26var.units = 'Sv (10^6 m^3/s)'
            # save variables
            Timevar[:] = Time[ncount-12:ncount]
            mocAtlantic26var[:] = mocAtlantic26[ncount-12:ncount]
            ncFile.close()
        else:
            ncount += 12
            # Read from file
            print '     Read previously computed Atlantic MOC at 26.5N'\
                  ' from file...'
            ncFile = netCDF4.Dataset(outputFileTseries, mode='r')
            Time[ncount-12:ncount] = ncFile.variables['Time'][:]
            mocAtlantic26[ncount-12:ncount] = \
                ncFile.variables['mocRapidArrayMax'][:]
            ncFile.close()

    # Store max Atlantic MOC at 26.5 to dictionary
    mocDictTseries = {'Time': {'dims': ('Time'), 'data': Time},
                      'mocRapidArrayMax': {'dims': ('Time'),
                                           'data': mocAtlantic26}}

    return (mocDictClimo, mocDictTseries)


# def _compute_moc_analysismember(config):
#
#     return (mocDictClimo, mocDictTseries)


def _compute_transport(maxEdgesInTransect, transectEdgeGlobalIDs,
                       transectEdgeMaskSigns, nz, dvEdge, refLayerThickness,
                       horizontalVel):

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
    return transportZ


def _compute_moc(latBins, nz, latCell, regionCellMask, transportZ, velArea):

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
    return mocTop

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
