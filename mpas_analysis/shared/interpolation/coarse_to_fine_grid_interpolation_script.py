#!/usr/bin/env python
"""
    Interpolation of initial conditions from coarse restart onto a finer grid.

    Currently implementation interpolates

        * temperature
        * salinity
        * layerThickness

    Phillip Wolfram, Mark Petersen
    01/20/2017
"""
import shutil
import numpy as np
import matplotlib.pyplot as plt
import netCDF4
from scipy.spatial import cKDTree as KDTree

# special timing function
# from http://stackoverflow.com/questions/5478351/python-time-measure-function
from contextlib import contextmanager
import time
@contextmanager
def timeit_context(name):
    startTime = time.time()
    yield
    elapsedTime = time.time() - startTime
    print('[{}] finished in {} s'.format(name, int(elapsedTime)))

def print_warning_xy(): #{{{
    print 'WARNING!!!:  Function currently uses xCell and yCell for '\
          'nearest neighbor interpolation.  latCell and lonCell may '\
          'better options.  This implementation is not strictly general '\
          'but may work in practice for a fully unstructured grid.'
    return #}}}

# general interpolation functions
def get_point_vectors3d(ncdfdata): #{{{
    """
    inputs:
        ncdfdata is an opened netCDF4.Dataset
    outputs:
        position tuple
    """
    maxLevelCell = ncdfdata.variables['maxLevelCell'][:]
    npoints = np.sum(maxLevelCell)

    # point data
    print_warning_xy()
    lonCell = ncdfdata.variables['xCell'][:]
    latCell = ncdfdata.variables['yCell'][:]
    # note this is a coarse proxy for the cell center, but it should work for
    # this application
    vertLevel = ncdfdata.variables['refBottomDepth'][:]

    x = np.zeros((npoints))
    y = np.zeros((npoints))
    z = np.zeros((npoints))

    n = 0
    for aCell, nLevel in enumerate(maxLevelCell):
        # vertical grid
        x[n:(n+nLevel)] = lonCell[aCell]
        y[n:(n+nLevel)] = latCell[aCell]
        z[n:(n+nLevel)] = vertLevel[:nLevel]
        # counter
        n += nLevel

    return np.vstack((x,y,z)).T #}}}

def get_point_vectors2d(ncdfdata): #{{{
    """
    inputs:
        ncdfdata is an opened netCDF4.Dataset
    outputs:
        position
    """
    # point data
    print_warning_xy()
    x = ncdfdata.variables['xCell'][:]
    y = ncdfdata.variables['yCell'][:]

    return np.vstack((x,y)).T #}}}

def get_3dcell_data(ncdfdata, datanames=[], ts=0): #{{{
    """
    inputs:
        ncdfdata is an opened netCDF4.Dataset
        datanames is a list of datanames to be interpolated
        ts is the time step number
    outputs:
        position and data
    """
    maxLevelCell = ncdfdata.variables['maxLevelCell'][:]
    npoints = np.sum(maxLevelCell)

    data = np.zeros((len(datanames),npoints))

    n = 0
    for aCell, nLevel in enumerate(maxLevelCell):
        # 3dcell data
        for ad, dataname in enumerate(datanames):
            data[ad,n:(n+nLevel)] = ncdfdata.variables[dataname][ts, aCell, :nLevel]
        # counter
        n += nLevel

    return data #}}}

def get_2dcell_data(ncdfdata, datanames=[], ts=0): #{{{
    """
    inputs:
        ncdfdata is an opened netCDF4.Dataset
        datanames is a list of datanames to be interpolated
        ts is the time step number
    outputs:
        position and data
    """
    nCells = len(ncdfdata.dimensions['nCells'])
    data = np.zeros((len(datanames),nCells))

    for ad, dataname in enumerate(datanames):
        if dataname == 'ssh':
            for ac, maxLevel in enumerate(ncdfdata.variables['maxLevelCell']):
                data[ad,ac] = -ncdfdata.variables['bottomDepth'][ac] + \
                        np.sum(ncdfdata.variables['layerThickness'][ts,ac,:maxLevel],axis=-1)


    return data #}}}

def set_3dcell_data(ncdfdata, datanames, data, ts=0): #{{{
    """
    inputs:
        ncdfdata is an opened netCDF4.Dataset
        datanames is a list of datanames
        data is datavector corresponding to datanames
        ts is the time step number
    outputs:
        in place=> updata data in ncdfdata for datanames
    """
    maxLevelCell = ncdfdata.variables['maxLevelCell'][:]
    vertLevel = ncdfdata.variables['refBottomDepth'][:]

    n = 0
    for aCell, nLevel in enumerate(maxLevelCell):
        for ad, dataname in enumerate(datanames):
            ncdfdata.variables[dataname][ts, aCell, :nLevel] = data[ad,n:(n+nLevel)]
        n+= nLevel

    return None #}}}

def grid_interp(cgrid, fgrid, ogrid, interiorscalars=['temperature','salinity']): #{{{
    """
    Perform interpolation from a coarse grid (cgrid) onto a find grid (fgrid),
    writing to an output file (ogrid). This interpolation interpolates scalars
    list from cgrid onto fgrid and uses nearest neighbor interpolation.

    Interpolated scalars are specified by interiorscalars and the
    layerThickness is interpolated on the output grid for z-star coordinates.

    Phillip Wolfram
    11/18/2015
    """

    with timeit_context('Make output file'):
        shutil.copyfile(fgrid, ogrid)

    with timeit_context('Open files'):
        cgriddata = netCDF4.Dataset(cgrid,'r')
        fgriddata = netCDF4.Dataset(fgrid,'r')
        ogriddata = netCDF4.Dataset(ogrid,'r+')

    with timeit_context('Build vectors of point locations x,y,z and data values'):
        cpos3d = get_point_vectors3d(cgriddata)
        cdata3dcell = get_3dcell_data(cgriddata, interiorscalars)
        fpos3d = get_point_vectors3d(fgriddata)

        cpos2d = get_point_vectors2d(cgriddata)
        cdata2dcell = get_2dcell_data(cgriddata, ['ssh'])
        fpos2d = get_point_vectors2d(fgriddata)

    with timeit_context('Building trees for coarse data'):
        tree3d = KDTree(cpos3d)
        tree2d = KDTree(cpos2d)

    with timeit_context('Finding nearest neighbors'):
        _, nearestneighs3d = tree3d.query(fpos3d)
        _, nearestneighs2d = tree2d.query(fpos2d)
    fdata3dcell = cdata3dcell[:,nearestneighs3d]
    fssh = cdata2dcell[:,nearestneighs2d]

    with timeit_context('Convert ssh to layerThickness'):
        rbd = fgriddata.variables['refBottomDepth'][:]
        relthick = np.hstack((rbd[0], np.diff(rbd)))
        layerThickness = relthick[np.newaxis,:] + (relthick/np.sum(relthick))[np.newaxis,:]*fssh.T

    with timeit_context('Saving data'):
        set_3dcell_data(ogriddata, interiorscalars, fdata3dcell)
        ogriddata.variables['layerThickness'][0,:,:] = layerThickness

    with timeit_context('Close files'):
        cgriddata.close()
        fgriddata.close()
        ogriddata.close()

    return None #}}}

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-c", "--coarsegrid", dest="cgrid", help="Input: Coarse grid to perform interpolation.")
    parser.add_argument("-f", "--finegrid", dest="fgrid", help="Input: Fine grid to perform interpolation.")
    parser.add_argument("-o", "--outputgrid", dest="ogrid", help="Output: Fine grid with interpolation of temperature and layerThickness.")

    args = parser.parse_args()

    if args.cgrid is None or args.fgrid is None or args.ogrid is None:
        parser.error('Must fully specify all inputs and outputs')
    if args.cgrid == args.ogrid:
        parser.error('Output grid cannot be specified as coarse input grid.')
    if args.fgrid == args.ogrid:
        parser.error('Output grid cannot be specified as fine input grid.')

    grid_interp(args.cgrid, args.fgrid, args.ogrid)
