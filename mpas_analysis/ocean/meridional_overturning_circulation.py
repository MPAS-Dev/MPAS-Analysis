import xarray as xr
import numpy as np
import datetime
import netCDF4

from ..shared.constants.constants import m3ps_to_Sv
from ..shared.plot.plotting import plot_moc

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    days_to_datetime

from ..shared.climatology import climatology

def moc_streamfunction(config, streamMap=None, variableMap=None): #{{{
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

    # read paramoeters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

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
    #outputTimes = config.getExpression(sectionName, 'comparisonTimes')

    #fileNameMask = config.get('moc_postprocess', 'filenamemask') # for regional
    #mocLat = config.getExpression('moc_postprocess', 'mocLat', usenumpy=True) # probably to be removed

    # open files, load variables
    ncFile = netCDF4.Dataset(restartFile, mode='r')
    dvEdge = ncFile.variables['dvEdge'][:]
    areaCell = ncFile.variables['areaCell'][:]
    refBottomDepth = ncFile.variables['refBottomDepth'][:]
    latCell = ncFile.variables['latCell'][:]
    ncFile.close()
    latCell = np.rad2deg(latCell)

    variableList = ['avgNormalVelocity',
                    'avgVertVelocityTop']
    ds = open_multifile_dataset(fileNames=inputFiles,
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=variableList,
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    # Compute annual climatology
    annualClimatology = climatology.compute_annual_climatology(ds, calendar)

    horizontalVel = annualClimatology.avgNormalVelocity
    verticalVel = annualClimatology.avgVertVelocityTop

    nVertLevels = annualClimatology.nVertLevels.size
    refLayerThickness = np.zeros(nVertLevels)
    refLayerThickness[0] = refBottomDepth[0]
    refLayerThickness[1:nVertLevels] = refBottomDepth[1:nVertLevels]-refBottomDepth[0:nVertLevels-1]

    refTopDepth = np.zeros(nVertLevels+1)
    refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]

    print '  Compute post-processed global MOC...'
    mocLat = np.arange(-90.0,90.1,1.) # 1deg meridional bins
    mocLat = mocLat*np.pi/180.0; # convert to radians
    velArea = verticalVel*areaCell[:,np.newaxis]
    #newlat = xr.DataArray(latCell, coords=[latCell], dims=['nCells'])
    #velArea.coords['latCell'] = ('nCells',newlat)
    #print 'velArea=\n',velArea
    #print ''
    #print 'velArea_binned.sum=\n',velArea.groupby_bins('latCell',mocLat).sum()
    #mocTop = velArea.groupby_bins('latCell',mocLat).sum()
    #print ''
    #print 'mocTop=\n',mocTop
    #print mocTop.cumsum('latCell_bin')
    nLat = np.size(mocLat)
    mocTop = np.zeros((nLat,nVertLevels+1))
    velArea = velArea.values
    for iLat in range(1,nLat):
        print "ilat=",iLat
    #    #indlat = np.logical_and(latCell >= mocLat[iLat-1], latCell <  mocLat[iLat])
        indlat = np.where( np.logical_and(latCell >= mocLat[iLat-1], latCell <  mocLat[iLat]) )
        indlat = np.squeeze(indlat)
        mocTop[iLat,:] = mocTop[iLat-1,:] + velArea[indlat,:].sum()
#        for k in range(1,nVertLevels+1):
#            mocTop[iLat,k] = mocTop[iLat-1,k] + verticalVel[indlat,k]*areaCell[indlat].sum()
#            print "k=",k
#            velArea = verticalVel[indlat,k]*areaCell[indlat]
#            #print velArea
#            mocTop[iLat,k] = mocTop[iLat-1,k] + velArea.sum('nCells')

        #for k in range(1,nVertLevels+1):
        #    #mocTop[iLat,k] = mocTop[iLat-1,k] + np.sum(verticalVel[indlat,k]*areaCell[indlat])
        #    mocTop[iLat,k] = mocTop[iLat-1,k] + verticalVel[indlat,k]*areaCell[indlat].sum()
    # convert m^3/s to Sverdrup
    mocTop = mocTop*m3ps_to_Sv
    print '  Plot post-processed global MOC...'
    title = 'Model global MOC'
    figureName = '{}/moc_global_{}_years{:04d}-{:04d}.png'.format(
                 plotsDirectory, mainRunName, startYear, endYear)
    plot_moc(config, mocLat, refTopDepth, mocTop, figureName, title)

#    print '  Compute post-processed regional MOC...'
#    regionalFile = '/global/project/projectdirs/acme/mapping/grids/EC60to30v1_SingleRegionAtlanticWTransportTransects_masks.nc'
#    ncFileRegional = netCDF4.Dataset(regionalFile, mode='r')
#    print '   Reading transect information...'
#    nTransects = ncFileRegional.dimensions['nTransects'].size
#    maxEdgesInTransect = ncFileRegional.dimensions['maxEdgesInTransect'].size
#    transectEdgeMaskSigns = ncFileRegional.variables['transectEdgeMaskSigns'][:,:]
#    transectEdgeGlobalIDs = ncFileRegional.variables['transectEdgeGlobalIDs'][:,:]
#    print '   Reading transect information...'
#    nRegions = ncFileRegional.dimensions['nRegions'].size
#    regionCellMasks = ncFileRegional.variables['regionCellMasks']
#    ncFileRegional.close()
#    print '   Compute volume transport through transects...'
#    transport = np.zeros(nTransects);
#    transportZ = np.zeros([nVertLevels,nTransects]);
#    transportZEdge = np.zeros([nVertLevels,maxEdgesInTransect,nTransects]);
#    for iTransect in range(nTransects):
#        for i in range(maxEdgesInTransect):
#            if transectEdgeGlobalIDs[iTransect, i]==0:
#                break
#            iEdge = transectEdgeGlobalIDs[iTransect, i] - 1 # subtract 1 because of python 0-indexing
#            for k in range(nVertLevels):
#                transportZEdge[k,i,iTransect] = ( transectEdgeMaskSigns[iEdge,iTransect]
#                        * horizontalVel[iEdge,k] * dvEdge[iEdge] * refLayerThickness[k]*m3ps_to_Sv )
#                transportZ[k,iTransect] = transportZ[k,iTransect] + transportZEdge[k,i,iTransect]
#                transport[iTransect] = transport[iTransect] + transportZEdge[k,i,iTransect]
#
#    mocLat = np.linspace(-34.5,70.0,210)
#    mocLat = mocLat*np.pi/180.0;
#    nLat = np.size(mocLat)
#    mocTop = np.zeros((nLat,nVertLevels+1))
#    print nRegions
#    for iRegion in range(nRegions):
#        # assume transects and regions have the same index ordering:
#        iTransect = iRegion
#
#        for k in range(1,nVertLevels+1):
#            mocTop[0,k] = mocTop[0,k-1] + transportZ[k-1,iTransect]/m3ps_to_Sv
#
#        for iLat in range(1,nLat):
#            ind = np.where( np.logical_and(np.logical_and(regionCellMasks[:,iRegion]==1,latCell >= mocLat[iLat-1] ),
#                                           latCell <  mocLat[iLat]) )
#            indlat = np.squeeze(indlat)
#
#            for k in range(1,nVertLevels+1):
#                mocTop[iLat,k] = mocTop[iLat-1,k] + np.sum(verticalVel[ind,k]*areaCell[ind])
#
#    # convert m^3/s to Sverdrup
#    mocTop = mocTop * m3ps_to_Sv
#
#    print '  Plot post-processed Atlantic MOC...'
#    title = 'Model Atlantic MOC'
#    figureName = '{}/moc_atlantic_{}_years{:04d}-{:04d}.png'.format(
#                 plotsDirectory, mainRunName, startYear, endYear)
#    plot_moc(config, mocLat, refTopDepth, mocTop, figureName, title)

    return #}}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
