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
03/28/2017
"""

import xarray as xr
import numpy as np
import netCDF4
import os
import warnings

from ..shared.constants.constants import m3ps_to_Sv, rad_to_deg
from ..shared.plot.plotting import plot_vertical_section,\
    timeseries_analysis_plot, setup_colormap

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    days_to_datetime

from ..shared.analysis_task import setup_task
from ..shared.climatology import climatology


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
    03/23/2017
    """

    # **** Initial settings ****

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    try:
        timeSeriesStatsMonthlyAnalysisMemberFlag = namelist.get(
            'config_am_timeseriesstatsmonthly_enable')
    except KeyError:
        warnings.warn('WARNING: namelist option config_am_timeseriesstatsmonthly_enable not found.\n'
                      'This likely indicates that the simulation you are analyzing was run with an\n'
                      'older version of MPAS-O that did not support this flag.  Assuming .true.')
        timeSeriesStatsMonthlyAnalysisMemberFlag = '.true.'

    if timeSeriesStatsMonthlyAnalysisMemberFlag == '.false.':
        raise RuntimeError('*** MPAS-Analysis relies on the existence of monthly\n'
                           '*** averaged files: make sure to enable the\n'
                           '*** timeSeriesStatsMonthly analysis member!')

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
    inputFilesTseries = historyStreams.readpath(streamName,
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
    try:
        mocAnalysisMemberFlag = namelist.get('config_am_mocstreamfunction_enable')
    except KeyError:
        warnings.warn('WARNING: namelist option config_am_mocstreamfunction_enable not found.\n'
                      'This likely indicates that the simulation you are analyzing was run with an\n'
                      'older version of MPAS-O that did not support this flag.  Assuming .false.')
        mocAnalysisMemberFlag = '.false.'

    # Check whether MOC Analysis Member is enabled
    if mocAnalysisMemberFlag == '.true.':
        # Add a moc_analisysMember_processing
        print '*** MOC Analysis Member is on ***'
        # (mocDictClimo, mocDictTseries) = _compute_moc_analysismember(config,
        #     streams, calendar, sectionName, dictClimo, dictTseries)
    else:
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
                 region, dictClimo['startYearClimo'], dictClimo['endYearClimo'],
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
    ds = xr.open_dataset(restartFile)
    latCell = ds.latCell
    latCell = latCell * rad_to_deg  # convert to degree
    refTopDepth = xr.DataArray(np.hstack((0, ds.refBottomDepth.values)),
                               coords=[ds.nVertLevelsP1])
    refLayerThickness = xr.DataArray(np.hstack((ds.refBottomDepth[0].values,
                                                ds.refBottomDepth.diff('nVertLevels').values)),
                                     coords=[ds.nVertLevels])

    return ds.dvEdge, ds.areaCell, ds.refBottomDepth, latCell, \
            ds. nVertLevels, ds.nVertLevelsP1, \
            refTopDepth, refLayerThickness  # }}}


def _compute_moc_climo_postprocess(config, runStreams, variableMap, calendar,
                                   sectionName, regionNames, dictClimo):  # {{{

    '''compute mean MOC streamfunction as a post-process'''

    dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        nVertLevelsP1, refTopDepth, refLayerThickness = _load_mesh(runStreams)

    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']

    # Load basin region related variables and save them to dictionary
    # NB: The following will need to change with new regional mapping files
    regionMaskFiles = config.get(sectionName, 'regionMaskFiles')
    if not os.path.exists(regionMaskFiles):
        raise IOError('Regional masking file for MOC calculation '
                      'does not exist')
    iRegion = 0
    dictRegion = {}
    for region in regionNames:
        print '\n  Reading region and transect mask for {}...'.format(region)
        dsRegion = xr.open_dataset(regionMaskFiles)
        transectEdgeMaskSigns = \
            dsRegion.transectEdgeMaskSigns.isel(nTransects=iRegion)
        regionCellMask = dsRegion.regionCellMasks.isel(nRegions=iRegion)
        iRegion += 1

        # this will only have the last entry in regionNames
        dictRegion.update({
            '{}CellMask'.format(region): regionCellMask,
            'transectEdgeMaskSigns{}'.format(region): transectEdgeMaskSigns})
    # Add Global regionCellMask=1 everywhere to make the algorithm
    # for the global moc similar to that of the regional moc
    dictRegion['GlobalCellMask'] = xr.ones_like(regionCellMask).rename('GlobalCellMask')
    regionNames[:0] = ['Global']

    # Compute and plot annual climatology of MOC streamfunction
    print '\n  Compute and/or plot post-processed MOC climatological '\
          'streamfunction...'
    simulationStartTime = get_simulation_start_time(runStreams)
    outputDirectory = build_config_full_path(config, 'output',
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
            config=config,
            simulationStartTime=simulationStartTime,
            timeVariableName='Time',
            variableList=variableList,
            variableMap=variableMap,
            startDate=dictClimo['startDateClimo'],
            endDate=dictClimo['endDateClimo'])

        changed, startYear, endYear = \
            climatology.update_start_end_year(ds, config, calendar)
        if changed:
            # update the file name in case the start and end years changed
            outputFileClimo = \
                '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
                    outputDirectory, startYear, endYear)

            dictClimo['startYearClimo'] = startYear
            dictClimo['endYearClimo'] = endYear

        # Compute annual climatology
        annualClimatology = ds.mean('Time')

        # Convert to numpy arrays
        horizontalVel = annualClimatology.avgNormalVelocity
        velArea = annualClimatology.avgVertVelocityTop * areaCell

        # Create dictionary for MOC climatology (NB: need this form
        # in order to convert it to xarray dataset later in the script)
        mocDictClimo = {'depth': {'dims': ('nz'), 'data': refTopDepth}}
        for region in regionNames:
            print '   Compute {} MOC...'.format(region)
            print '    Compute transport through region southern transect...'
            if region == 'Global':
                transportZ = xr.zeros_like(nVertLevels).rename('transportZ')
            else:
                transectEdgeMaskSigns = \
                    dictRegion['transectEdgeMaskSigns{}'.format(region)]
                transportZ = _compute_transport(transectEdgeMaskSigns,
                                                nVertLevels, dvEdge,
                                                refLayerThickness,
                                                horizontalVel)

            regionCellMask = dictRegion['{}CellMask'.format(region)]
            latBinSize = config.getExpression(sectionName,
                                              'latBinSize{}'.format(region))
            if region == 'Global':
                latBins = np.arange(-90.0, 90.1, latBinSize)
            else:
                cellMask = dictRegion['{}CellMask'.format(region)]
                latBins = latCell.where(cellMask)
                latBins = np.arange(np.amin(latBins),
                                    np.amax(latBins)+latBinSize,
                                    latBinSize)
            latBins = xr.DataArray(latBins, coords=[latBins], dims=['latBins'])
            mocTop = _compute_moc(latBins, nVertLevelsP1, latCell,
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

    def write_file(dsMOCTimeSeries):
        # make sure the Time dimension is in ascending order by Time
        indices = dsMOCTimeSeries.Time.argsort()
        dsMOCTimeSeries = dsMOCTimeSeries.isel(Time=indices)
        dsMOCTimeSeries.Time.attrs['units'] = 'days since 0001-01-01'
        dsMOCTimeSeries.mocAtlantic26.attrs['units'] = 'Sv (10^6 m^3/s)'
        dsMOCTimeSeries.mocAtlantic26.attrs['description'] = \
            ('Max MOC Atlantic streamfunction nearest to '
             'RAPID Array latitude (26.5N)')
        dsMOCTimeSeries.to_netcdf(outputFileTseries)
        return dsMOCTimeSeries

    # Compute and plot time series of Atlantic MOC at 26.5N (RAPID array)
    print '\n  Compute and/or plot post-processed Atlantic MOC '\
          'time series...'
    print '   Load data...'

    simulationStartTime = get_simulation_start_time(runStreams)
    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']

    dvEdge, areaCell, refBottomDepth, latCell, nVertLevels, \
        nVertLevelsP1, refTopDepth, refLayerThickness = _load_mesh(runStreams)

    ds = open_multifile_dataset(fileNames=dictTseries['inputFilesTseries'],
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=dictTseries['startDateTseries'],
                                endDate=dictTseries['endDateTseries'])
    latAtlantic = mocDictClimo['latAtlantic']['data']

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
        # read in what we have so far
        dsMOCTimeSeries = xr.open_dataset(outputFileTseries,
                                          decode_times=False)
        # force loading and then close so we can overwrite the file later
        dsMOCTimeSeries.load()
        dsMOCTimeSeries.close()
    else:
        # create a new data set from a dictionary
        dictonary = {
            'Time': {'dims': ('Time'),
                     'data': []},
            'mocAtlantic26': {'dims': ('Time'),
                              'data': []}}
        dsMOCTimeSeries = xr.Dataset.from_dict(dictonary)

    newEntries = False
    first = True
    for tIndex in range(ds.Time.size):
        if ds.Time[tIndex].values in dsMOCTimeSeries.Time.values:
            # we have it already, so no need to recompute
            continue
        if first:
            print '   Process and save MOC time series'
            first = False
        date = days_to_datetime(ds.Time[tIndex], calendar=calendar)
        print '     date: {:04d}-{:02d}'.format(date.year, date.month)
        horizontalVel = ds.avgNormalVelocity.isel(Time=tIndex)
        velArea = ds.avgVertVelocityTop.isel(Time=tIndex) * areaCell
        transportZ = _compute_transport(transectEdgeMaskSigns,
                                        nVertLevels, dvEdge,
                                        refLayerThickness,
                                        horizontalVel)
        mocTop = _compute_moc(latAtlantic, nVertLevelsP1, latCell,
                              regionCellMask, transportZ, velArea)
        mocAtlantic26 = mocTop.sel(latBins=26.5, method='nearest').max()

        dictonary = {'Time': {'dims': ('Time'), 'data': [ds.Time[tIndex]]},
                     'mocAtlantic26': {'dims': ('Time'),
                                       'data': [mocAtlantic26]}}
        # add the new time entry to the data set
        dsMOCTimeSeries = xr.concat([dsMOCTimeSeries,
                                     xr.Dataset.from_dict(dictonary)],
                                    dim='Time')
        newEntries = True
        if date.month == 12:
            # each year, save the data set for safe keeping
            print '   saving progress.'
            dsMOCTimeSeries = write_file(dsMOCTimeSeries)
            newEntries = False

    if newEntries:
        print '   saving final result.'
        dsMOCTimeSeries = write_file(dsMOCTimeSeries)

    return dsMOCTimeSeries  # }}}


# def _compute_moc_analysismember(config):
#
#     return (mocDictClimo, mocDictTseries)


def _compute_transport(transectEdgeMaskSigns, nz, dvEdge, refLayerThickness,
                       horizontalVel):  # {{{

    '''compute mass transport across southern transect of ocean basin'''

    transect = transectEdgeMaskSigns != 0
    transportZ = (refLayerThickness*horizontalVel.where(transect)
                  *transectEdgeMaskSigns.where(transect)
                  *dvEdge.where(transect)
                  ).sum('nEdges')

    return transportZ  # }}}


def _compute_moc(latBins, nz, latCell, regionCellMask, transportZ,
                 velArea):  # {{{

    '''compute meridionally integrated MOC streamfunction'''
    # note computation should be more cleanly facilitated 
    # via a multidimensional groupby_bins when available
    # (see https://github.com/pydata/xarray/pull/924)

    mocTop = xr.DataArray(np.zeros([len(nz), len(latBins)]),
                          coords=[nz, latBins])
    mocTop[1:, 0] = transportZ.cumsum(axis=0)
    for iLat in np.arange(1, len(latBins)):
        mask = np.logical_and(
                np.logical_and(latCell >= latBins[iLat-1],
                               latCell < latBins[iLat]),
                               regionCellMask > 0)
        mocTop[:, iLat] = velArea.where(mask).sum('nCells')
    mocTop = mocTop.cumsum('latBins')
    # convert m^3/s to Sverdrup
    mocTop *= m3ps_to_Sv
    return mocTop  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
