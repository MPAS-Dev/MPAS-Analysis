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
03/13/2017
"""

import matplotlib.pyplot as plt
import matplotlib.colors as cols

import xarray as xr
import datetime
import numpy as np
import netCDF4
import os

from ..shared.constants.constants import m3ps_to_Sv
from ..shared.plot.plotting import plot_vertical_section

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
    post-processing if not. Plots streamfunction as well as time series
    of max MOC at 24N.

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Authors: Mark Petersen, Phillip J. Wolfram, Milena Veneziani
    Modified: 02/23/2017
    """

    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    timeSeriesStatsMonthlyAnalysisMemberFlag = namelist.get(
        'config_am_timeseriesstatsmonthly_enable')
    if timeSeriesStatsMonthlyAnalysisMemberFlag == '.false.':
        raise RuntimeError('*** MPAS-Analysis relies on the existence of '
                           'monthly\n*** averaged files: make sure to '
                           'enable the\n*** timeSeriesStatsMonthly '
                           'analysis member!')

    mocAnalysisMemberFlag = namelist.get(
        'config_am_mocstreamfunction_enable')

    calendar = namelist.get('config_calendar_type')
    simulationStartTime = get_simulation_start_time(streams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    inputFiles = streams.readpath(streamName, startDate=startDate,
                                  endDate=endDate, calendar=calendar)
    print 'Reading files {} through {}'.format(inputFiles[0], inputFiles[-1])

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')

    try:
        restartFile = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for MOC calculation')

    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')

    sectionName = 'streamfunctionMOC'

    # Plotting variables
    colormapName = plt.get_cmap(config.get(sectionName, 'colormapName'))
    colormapIndices = config.getExpression(sectionName, 'colormapIndices')
    colormapName = cols.ListedColormap(colormapName(colormapIndices),
                                       'colormapName')
    colorbarLabel = '[Sv]'

    # Check whether MOC Analysis Member is enabled
    if mocAnalysisMemberFlag == '.true.':

        # Add a moc_analisysMember_processing
        print '*** MOC Analysis Member is on ***'

    else:

        outputDirectory = buildConfigFullPath(config, 'output',
                                              'mpasClimatologySubdirectory')
        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass
        outputFile = '{}/mocStreamfunction_years{:04d}-{:04d}.nc'.format(
                      outputDirectory, startYear, endYear)

        if not os.path.exists(outputFile):
            # Load mesh related variables
            ncFile = netCDF4.Dataset(restartFile, mode='r')
            dvEdge = ncFile.variables['dvEdge'][:]
            areaCell = ncFile.variables['areaCell'][:]
            refBottomDepth = ncFile.variables['refBottomDepth'][:]
            latCell = ncFile.variables['latCell'][:]
            latCell = latCell * 180.0 / np.pi  # convert to degree
            ncFile.close()
            nVertLevels = len(refBottomDepth)
            refTopDepth = np.zeros(nVertLevels+1)
            refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]

            refLayerThickness = np.zeros(nVertLevels)
            refLayerThickness[0] = refBottomDepth[0]
            refLayerThickness[1:nVertLevels] =  \
                refBottomDepth[1:nVertLevels] - \
                refBottomDepth[0:nVertLevels-1]

            variableList = ['avgNormalVelocity',
                            'avgVertVelocityTop']
            ds = open_multifile_dataset(
                fileNames=inputFiles,
                calendar=calendar,
                simulationStartTime=simulationStartTime,
                timeVariableName='Time',
                variableList=variableList,
                variableMap=variableMap,
                startDate=startDate,
                endDate=endDate)

            # Compute annual climatology
            annualClimatology = climatology.compute_annual_climatology(
                                ds, calendar)

            horizontalVel = annualClimatology.avgNormalVelocity
            verticalVel = annualClimatology.avgVertVelocityTop
            # Convert to numpy
            horizontalVel = horizontalVel.values
            verticalVel = verticalVel.values
            velArea = verticalVel * areaCell[:, np.newaxis]

            print '\n  Compute post-processed global MOC...'
            mocLatGlobal = np.arange(-90.0, 90.1, 1.)  # 1deg meridional bins
            nLat = np.size(mocLatGlobal)
            mocTop = np.zeros((nLat, nVertLevels+1))
            for iLat in range(1, nLat):
                indlat = np.logical_and(latCell >= mocLatGlobal[iLat-1],
                                        latCell < mocLatGlobal[iLat])
                mocTop[iLat, :] = mocTop[iLat-1, :] + \
                                  velArea[indlat, :].sum(axis=0)
            # convert m^3/s to Sverdrup
            mocGlobal = mocTop * m3ps_to_Sv
            mocGlobal = mocGlobal.T

            print '\n  Compute post-processed regional MOC...'
            # The following will need to change with new regional mapping files
            regionMaskFiles = config.get(sectionName, 'regionMaskFiles')
            if not os.path.exists(regionMaskFiles):
                raise IOError('Regional masking file for MOC calculation '
                              'does not exist')
            #for iRegion in range(len(regionMaskFiles)):
            print '   Reading region and transect mask...'
            ncFileRegional = netCDF4.Dataset(regionMaskFiles, mode='r')
            iTransect = 0  # Atlantic
            maxEdgesInTransect = ncFileRegional.dimensions['maxEdgesInTransect'].size
            transectEdgeMaskSigns = ncFileRegional.variables['transectEdgeMaskSigns'][:, iTransect]
            transectEdgeGlobalIDs = ncFileRegional.variables['transectEdgeGlobalIDs'][iTransect, :]
            #nRegions = ncFileRegional.dimensions['nRegions'].size
            regionCellMask = ncFileRegional.variables['regionCellMasks'][:, 0]
            ncFileRegional.close()

            print '   Compute volume transport through southern transect...'
            transportZEdge = np.zeros([nVertLevels, maxEdgesInTransect])
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

            print '   Compute regional MOC...'
            indRegion = np.where(regionCellMask == 1)
            mocLatRegional = latCell[indRegion]
            mocLatRegional = np.arange(np.amin(mocLatRegional)-0.5,
                                       np.amax(mocLatRegional)+0.5, 0.5)
            nLat = np.size(mocLatRegional)
            mocTop = np.zeros((nLat, nVertLevels+1))
            mocTop[0, range(1, nVertLevels+1)] = transportZ.cumsum()
            for iLat in range(1, nLat):
                indlat = np.logical_and(np.logical_and(
                       regionCellMask == 1, latCell >= mocLatRegional[iLat-1]),
                       latCell < mocLatRegional[iLat])
                mocTop[iLat, :] = mocTop[iLat-1, :] + \
                                  velArea[indlat, :].sum(axis=0)
            # convert m^3/s to Sverdrup
            mocRegional = mocTop * m3ps_to_Sv
            mocRegional = mocRegional.T

            # Save to file
            ncFile = netCDF4.Dataset(outputFile, mode='w')
            # create dimensions
            ncFile.createDimension('nxGlobal', len(mocLatGlobal))
            ncFile.createDimension('nxAtlantic', len(mocLatRegional))
            ncFile.createDimension('nz', len(refTopDepth))
            # create variables
            latBinGlobal = ncFile.createVariable(
                'latBinGlobal', 'f4', ('nxGlobal',))
            latBinGlobal.description = \
                'latitude bins for MOC Global streamfunction'
            latBinGlobal.units = 'degrees (-90 to 90)'
            latBinAtlantic = ncFile.createVariable(
                'latBinAtlantic', 'f4', ('nxAtlantic',))
            latBinAtlantic.description = \
                'latitude bins for MOC Atlantic streamfunction'
            latBinAtlantic.units = 'degrees (-90 to 90)'
            depth = ncFile.createVariable('depth', 'f4', ('nz',))
            depth.description = 'depth'
            depth.units = 'meters'
            fldGlobal = ncFile.createVariable(
                'mocGlobal', 'f4', ('nz', 'nxGlobal'))
                #'mocGlobal', 'f4', ('nxGlobal', 'nz'))
            fldGlobal.description = \
                'MOC Global streamfunction, annual climatology'
            fldGlobal.units = 'Sv (10^6 m^3/s)'
            fldAtlantic = ncFile.createVariable(
                'mocAtlantic', 'f4', ('nz', 'nxAtlantic'))
                #'mocAtlantic', 'f4', ('nxAtlantic', 'nz'))
            fldAtlantic.description = \
                'MOC Atlantic streamfunction, annual climatology'
            fldAtlantic.units = 'Sv (10^6 m^3/s)'
            # save variables
            latBinGlobal[:] = mocLatGlobal
            latBinAtlantic[:] = mocLatRegional
            depth[:] = refTopDepth
            fldGlobal[:, :] = mocGlobal
            fldAtlantic[:, :] = mocRegional
            ncFile.close()
        else:
            # Read from file
            ncFile = netCDF4.Dataset(outputFile, mode='r')
            mocLatGlobal = ncFile.variables['latBinGlobal'][:]
            mocLatRegional = ncFile.variables['latBinAtlantic'][:]
            refTopDepth = ncFile.variables['depth'][:]
            mocGlobal = ncFile.variables['mocGlobal'][:, :]
            mocRegional = ncFile.variables['mocAtlantic'][:, :]
            ncFile.close()

        print '\n  Plot post-processed global MOC...'
        title = 'Global MOC (ANN, years {:04d}-{:04d})\n {}'.format(
                 startYear, endYear, mainRunName)
        regionName = 'global'
        figureName = '{}/moc_{}_{}_years{:04d}-{:04d}.png'.format(
                      plotsDirectory, regionName, mainRunName,
                      startYear, endYear)
        contourLevels = config.getExpression(
                            sectionName, '{}ContourLevels'.format(regionName),
                            usenumpyfunc=True)
        colorbarLevels = config.getExpression(
                            sectionName, '{}ColorbarLevels'.format(regionName))

        plot_vertical_section(config, mocLatGlobal, refTopDepth, mocGlobal,
                              colormapName, colorbarLevels, contourLevels,
                              colorbarLabel, title, 'latitude [deg]',
                              figureName)

        print '  Plot post-processed Atlantic MOC...'
        title = 'Atlantic MOC (ANN, years {:04d}-{:04d})\n {}'.format(
                 startYear, endYear, mainRunName)
        regionName = 'atlantic'
        figureName = '{}/moc_{}_{}_years{:04d}-{:04d}.png'.format(
                      plotsDirectory, regionName, mainRunName,
                      startYear, endYear)
        contourLevels = config.getExpression(
                           sectionName, '{}ContourLevels'.format(regionName),
                           usenumpyfunc=True)
        colorbarLevels = config.getExpression(
                           sectionName, '{}ColorbarLevels'.format(regionName))

        plot_vertical_section(config, mocLatRegional, refTopDepth, mocRegional,
                              colormapName, colorbarLevels, contourLevels,
                              colorbarLabel, title, 'latitude [deg]',
                              figureName)

    return  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
