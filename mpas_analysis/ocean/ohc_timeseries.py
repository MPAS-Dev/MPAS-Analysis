import numpy as np
import netCDF4
from netCDF4 import Dataset as netcdf_dataset
import xarray as xr
import pandas as pd
import datetime

from ..shared.mpas_xarray.mpas_xarray import preprocess_mpas, \
    remove_repeated_time_index

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.timekeeping.Date import Date


def ohc_timeseries(config, streamMap=None, variableMap=None):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.
    config is an instance of an MpasAnalysisConfigParser containing
    configuration options.

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani
    Last Modified: 02/02/2017
    """

    # read parameters from config file
    mainRunName = config.get('runs', 'mainRunName')
    preprocessedReferenceRunName = config.get('runs',
                                              'preprocessedReferenceRunName')
    preprocessedInputDirectory = config.get('oceanPreprocessedReference',
                                            'baseDirectory')

    compareWithObservations = config.getboolean('timeSeriesOHC',
                                                'compareWithObservations')

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')

    yearOffset = config.getint('time', 'yearOffset')

    movingAveragePoints = config.getint('timeSeriesOHC', 'movingAveragePoints')

    regions = config.getExpression('regions', 'regions')
    plotTitles = config.getExpression('regions', 'plotTitles')
    regionIndicesToPlot = config.getExpression('timeSeriesOHC',
                                               'regionIndicesToPlot')

    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    # Note: input file, not a mesh file because we need dycore specific fields
    # such as refBottomDepth and namelist fields such as config_density0, as
    # well as simulationStartTime, that are not guaranteed to be in the mesh
    # file.
    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for OHC calculation')

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    inFiles = streams.readpath(streamName, startDate=startDate,
                               endDate=endDate)
    print 'Reading files {} through {}'.format(inFiles[0], inFiles[-1])

    # Define/read in general variables
    print '  Read in depth and compute specific depth indexes...'
    ncFile = netcdf_dataset(restartFile, mode='r')
    # reference depth [m]
    depth = ncFile.variables['refBottomDepth'][:]
    # simulation start time
    simulationStartTime = netCDF4.chartostring(
        ncFile.variables['simulationStartTime'][:])
    simulationStartTime = str(simulationStartTime)
    ncFile.close()
    # specific heat [J/(kg*degC)]
    cp = namelist.getfloat('config_specific_heat_sea_water')
    # [kg/m3]
    rho = namelist.getfloat('config_density0')
    factor = 1e-22*rho*cp

    k700m = np.where(depth > 700.)[0][0] - 1
    k2000m = np.where(depth > 2000.)[0][0] - 1

    kbtm = len(depth)-1

    # Load data
    print '  Load ocean data...'
    variableList = ['avgLayerTemperature',
                    'sumLayerMaskValue',
                    'avgLayerArea',
                    'avgLayerThickness']
    ds = xr.open_mfdataset(
        inFiles,
        preprocess=lambda x: preprocess_mpas(x,
                                             yearoffset=yearOffset,
                                             timestr='Time',
                                             onlyvars=variableList,
                                             varmap=variableMap))

    ds = remove_repeated_time_index(ds)

    # convert the start and end dates to datetime objects using
    # the Date class, which ensures the results are within the
    # supported range
    timeStart = Date(startDate).to_datetime(yearOffset)
    timeEnd = Date(endDate).to_datetime(yearOffset)
    # select only the data in the specified range of years
    ds = ds.sel(Time=slice(timeStart, timeEnd))

    # Select year-1 data and average it (for later computing anomalies)
    timeStartFirstYear = Date(simulationStartTime).to_datetime(yearOffset)
    if timeStartFirstYear < timeStart:
        startDateFirstYear = simulationStartTime
        endDateFirstYear = '{}-12-31_23:59:59'.format(startDateFirstYear[0:4])
        filesFirstYear = streams.readpath(streamName,
                                          startDate=startDateFirstYear,
                                          endDate=endDateFirstYear)
        dsFirstYear = xr.open_mfdataset(
            filesFirstYear,
            preprocess=lambda x: preprocess_mpas(x,
                                                 yearoffset=yearOffset,
                                                 timestr='Time',
                                                 onlyvars=variableList,
                                                 varmap=variableMap))

        dsFirstYear = remove_repeated_time_index(dsFirstYear)
    else:
        timeStart = datetime.datetime(timeStart.year, 1, 1)
        timeEnd = datetime.datetime(timeStart.year, 12, 31)
        dsFirstYear = ds.sel(Time=slice(timeStart, timeEnd))
    meanFirstYear = dsFirstYear.mean('Time')

    print '  Compute temperature anomalies...'
    avgLayerTemperature = ds.avgLayerTemperature
    avgLayerTemperatureFirstYear = meanFirstYear.avgLayerTemperature

    avgLayTemperatureAnomaly = (avgLayerTemperature -
                                avgLayerTemperatureFirstYear)

    yearStart = (pd.to_datetime(ds.Time.min().values)).year
    yearEnd = (pd.to_datetime(ds.Time.max().values)).year
    timeStart = datetime.datetime(yearStart, 1, 1)
    timeEnd = datetime.datetime(yearEnd, 12, 31)

    if preprocessedReferenceRunName != 'None':
        print '  Load in OHC from preprocessed reference run...'
        inFilesPreprocesses = '{}/OHC.{}.year*.nc'.format(
            preprocessedInputDirectory, preprocessedReferenceRunName)
        dsPreprocessed = xr.open_mfdataset(
            inFilesPreprocesses,
            preprocess=lambda x: preprocess_mpas(x, yearoffset=yearOffset))
        dsPreprocessed = remove_repeated_time_index(dsPreprocessed)
        yearEndPreprocessed = \
            (pd.to_datetime(dsPreprocessed.Time.max().values)).year
        if yearStart <= yearEndPreprocessed:
            dsPreprocessedTimeSlice = dsPreprocessed.sel(Time=slice(timeStart,
                                                                    timeEnd))
        else:
            print '   Warning: Preprocessed time series ends before the ' \
                'timeSeries startYear and will not be plotted.'
            preprocessedReferenceRunName = 'None'

    sumLayerMaskValue = ds.sumLayerMaskValue
    avgLayerArea = ds.avgLayerArea
    avgLayerThickness = ds.avgLayerThickness

    print '  Compute OHC and make plots...'
    for index in range(len(regionIndicesToPlot)):
        regionIndex = regionIndicesToPlot[index]

        # Compute volume of each layer in the region:
        layerArea = sumLayerMaskValue[:, regionIndex, :] * \
            avgLayerArea[:, regionIndex, :]
        layerVolume = layerArea * avgLayerThickness[:, regionIndex, :]

        # Compute OHC:
        ohc = layerVolume * avgLayTemperatureAnomaly[:, regionIndex, :]
        # OHC over 0-bottom depth range:
        ohcTotal = ohc.sum('nVertLevels')
        ohcTotal = factor*ohcTotal

        # OHC over 0-700m depth range:
        ohc700m = factor*ohc[:, 0:k700m].sum('nVertLevels')

        # OHC over 700m-2000m depth range:
        ohc2000m = factor*ohc[:, k700m+1:k2000m].sum('nVertLevels')

        # OHC over 2000m-bottom depth range:
        ohcBottom = ohc[:, k2000m+1:kbtm].sum('nVertLevels')
        ohcBottom = factor*ohcBottom

        title = 'OHC, {}, 0-bottom (thick-), 0-700m (thin-), 700-2000m (--),' \
                ' 2000m-bottom (-.) \n {}'.format(plotTitles[regionIndex],
                                                  mainRunName)

        xLabel = 'Time [years]'
        yLabel = '[x$10^{22}$ J]'

        if preprocessedReferenceRunName != 'None':
            figureName = '{}/ohc_{}_{}_{}.png'.format(
                    plotsDirectory,
                    regions[regionIndex],
                    mainRunName,
                    preprocessedReferenceRunName)
            ohcPreprocessedTotal = dsPreprocessedTimeSlice.ohc_tot
            ohcPreprocessed700m = dsPreprocessedTimeSlice.ohc_700m
            ohcPreprocessed2000m = dsPreprocessedTimeSlice.ohc_2000m
            ohcPreprocessedBottom = dsPreprocessedTimeSlice.ohc_btm
            title = '{} (r), {} (b)'.format(title,
                                            preprocessedReferenceRunName)
            timeseries_analysis_plot(config, [ohcTotal, ohc700m, ohc2000m,
                                              ohcBottom, ohcPreprocessedTotal,
                                              ohcPreprocessed700m,
                                              ohcPreprocessed2000m,
                                              ohcPreprocessedBottom],
                                     movingAveragePoints, title,
                                     xLabel, yLabel, figureName,
                                     lineStyles=['r-', 'r-', 'r--', 'r-.',
                                                 'b-', 'b-', 'b--', 'b-.'],
                                     lineWidths=[2, 1, 1.5, 1.5, 2, 1, 1.5,
                                                 1.5])

        if (not compareWithObservations and
                preprocessedReferenceRunName == 'None'):
            figureName = '{}/ohc_{}_{}.png'.format(plotsDirectory,
                                                   regions[regionIndex],
                                                   mainRunName)
            timeseries_analysis_plot(config, [ohcTotal, ohc700m, ohc2000m,
                                              ohcBottom],
                                     movingAveragePoints, title,
                                     xLabel, yLabel, figureName,
                                     lineStyles=['r-', 'r-', 'r--', 'r-.'],
                                     lineWidths=[2, 1, 1.5, 1.5])
