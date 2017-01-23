#!/usr/bin/env python
"""
Compute and plot moc, as a post-processing computation

Author: Mark Petersen, Phillip J. Wolfram
Last Modified: 01/20/2017
"""

from netCDF4 import Dataset
import numpy as np
from ..shared.constants.constants import m3ps_to_Sv
from ..shared.plot.plotting import plot_moc

def moc_postprocess(config, streamMap=None, variableMap=None): #{{{
    """
    Plots a comparison of ACME/MPAS output to SST or MLD observations

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Authors: Mark Petersen, Phillip J. Wolfram
    Modified: 01/20/2017
    """

    # read parameters from config file
    fileName = config.get('moc_postprocess', 'filename')
    fileNameMesh = config.get('moc_postprocess', 'filenamemesh')
    fileNameMask = config.get('moc_postprocess', 'filenamemask')
    timeIndex = config.get('moc_postprocess', 'timeIndex')
    horVelName = config.get('moc_postprocess', 'horVelName')
    vertVelName = config.get('moc_postprocess', 'vertVelname')
    mocLat = config.getExpression('moc_postprocess', 'mocLat', usenumpy=True)

    # open files, load variables
    print('** load variables: \n')

    print 'Reading data from: ' + fileName
    ncfile = Dataset(fileName,'r')
    nCells = ncfile.dimensions['nCells'].size
    nVertLevels = ncfile.dimensions['nVertLevels'].size

    horVel = ncfile.variables[horVelName][:,:,:]
    vertVel = ncfile.variables[vertVelName][:,:,:]

    print 'Reading data from: ' + fileNameMesh
    ncfileMesh = Dataset(fileNameMesh,'r')
    dvEdge = ncfileMesh.variables['dvEdge'][:]
    areaCell = ncfileMesh.variables['areaCell'][:]
    latCell = ncfileMesh.variables['latCell'][:]*180./np.pi
    refBottomDepth = ncfileMesh.variables['refBottomDepth'][:]

    refLayerThickness = np.zeros(nVertLevels)
    refLayerThickness[0] = refBottomDepth[0]
    refLayerThickness[1:nVertLevels] = refBottomDepth[1:nVertLevels]-refBottomDepth[0:nVertLevels-1]

    refTopDepth = np.zeros(nVertLevels+1)
    refTopDepth[1:nVertLevels+1] = refBottomDepth[0:nVertLevels]

    print 'Reading transects from: ' + fileNameMask
    ncfileMask = Dataset(fileNameMask,'r')
    nTransects = ncfileMask.dimensions['nTransects'].size
    maxEdgesInTransect = ncfileMask.dimensions['maxEdgesInTransect'].size

    transectEdgeMaskSigns = ncfileMask.variables['transectEdgeMaskSigns'][:,:]
    transectEdgeGlobalIDs = ncfileMask.variables['transectEdgeGlobalIDs'][:,:]

    print 'Reading regionss from: ' + fileNameMask
    ncfileMask = Dataset(fileNameMask,'r')
    nRegions = ncfileMask.dimensions['nRegions'].size
    regionCellMasks = ncfileMask.variables['regionCellMasks']

    # compute transport through transects
    print('** compute transport: \n')

    # the volume transport
    transport = np.zeros(nTransects);
    transportZ = np.zeros([nVertLevels,nTransects]);
    transportZEdge = np.zeros([nVertLevels,maxEdgesInTransect,nTransects]);

    for iTransect in range(nTransects):
        for i in range(maxEdgesInTransect):
            if transectEdgeGlobalIDs[iTransect, i]==0:
                break
            iEdge = transectEdgeGlobalIDs[iTransect, i] - 1 # subtract 1 because of python 0-indexing
            for k in range(nVertLevels):
                transportZEdge[k,i,iTransect] = ( transectEdgeMaskSigns[iEdge,iTransect]
                        * horVel[timeIndex,iEdge,k] * dvEdge[iEdge] * refLayerThickness[k]*m3ps_to_Sv )
                transportZ[k,iTransect] = transportZ[k,iTransect] + transportZEdge[k,i,iTransect]
                transport[iTransect] = transport[iTransect] + transportZEdge[k,i,iTransect]

    # compute MOC
    print('** compute moc: \n')

    nLat = np.size(mocLat)
    mocTop = np.zeros((nLat,nVertLevels+1))

    for iRegion in range(nRegions):
        # assume transects and regions have the same index ordering:
       iTransect = iRegion

       for k in range(1,nVertLevels+1):
           mocTop[0,k] = mocTop[0,k-1] + transportZ[k-1,iTransect]/m3ps_to_Sv

       for iLat in range(1,nLat):
           ind =np.logical_and(np.logical_and(regionCellMasks[:,iRegion]==1,latCell >= mocLat[iLat-1] ), latCell <  mocLat[iLat])

           for k in range(1,nVertLevels+1):
               mocTop[iLat,k] = mocTop[iLat-1,k] + np.sum(vertVel[timeIndex,ind,k]*areaCell[ind])

       # convert m^3/s to Sverdrup
       mocTop = mocTop * m3ps_to_Sv

    # plot MOC
    print('** plot moc: \n')
    plot_moc(config, mocLat, refTopDepth, mocTop)

    ncfile.close()
    ncfileMask.close()

    return #}}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
