'''
Functions for creating SCRIP files, used to create mapping files

Functions
---------
mpas_file_to_scrip - create a SCRIP file for an MPAS mesh

lat_lon_file_to_scrip - create a SCRIP file based on a lat-lon grid from a file

lat_lon_array_to_scrip - create a SCRIP file based on lat and lon arrays

Author
------
Xylar Asay-Davis

Last Modified
-------------
02/23/2017
'''

import netCDF4
import numpy
import sys


def mpas_file_to_scrip(mpasFileName, scripFileName):  # {{{
    '''
    Given an MPAS mesh file, create a SCRIP file based on the mesh.

    Parameters
    ----------
    mpasFileName : str
        The path of the file containing the source MPAS mesh

    scripFileName : str
        The path to which the SCRIP file should be written

    Authors
    ------
    Xylar Asay-Davis, Matthew Hoffman, Douglas Jacobsen

    Last Modified
    -------------
    02/20/2017
    '''
    inFile = netCDF4.Dataset(mpasFileName, 'r')
    outFile = netCDF4.Dataset(scripFileName, 'w')

    # Get info from input file
    latCell = inFile.variables['latCell'][:]
    lonCell = inFile.variables['lonCell'][:]
    latVertex = inFile.variables['latVertex'][:]
    lonVertex = inFile.variables['lonVertex'][:]
    verticesOnCell = inFile.variables['verticesOnCell'][:]
    nEdgesOnCell = inFile.variables['nEdgesOnCell'][:]
    nCells = len(inFile.dimensions['nCells'])
    maxVertices = len(inFile.dimensions['maxEdges'])
    areaCell = inFile.variables['areaCell'][:]
    sphereRadius = float(inFile.sphere_radius)

    _create_scrip(outFile, grid_size=nCells, grid_corners=maxVertices,
                  grid_rank=1, units='radians')

    grid_area = outFile.createVariable('grid_area', 'f8', ('grid_size',))
    grid_area.units = 'radian^2'
    # SCRIP uses square radians
    grid_area[:] = areaCell[:] / (sphereRadius**2)

    outFile.variables['grid_center_lat'][:] = latCell[:]
    outFile.variables['grid_center_lon'][:] = lonCell[:]
    outFile.variables['grid_dims'][:] = nCells
    outFile.variables['grid_imask'][:] = 1

    # grid corners:
    grid_corner_lon = numpy.zeros((nCells, maxVertices))
    grid_corner_lat = numpy.zeros((nCells, maxVertices))
    for iVertex in range(maxVertices):
        cellIndices = numpy.arange(nCells)
        # repeat the last vertex wherever iVertex > nEdgesOnCell
        localVertexIndices = numpy.minimum(nEdgesOnCell-1, iVertex)
        vertexIndices = verticesOnCell[cellIndices, localVertexIndices] - 1
        grid_corner_lat[cellIndices, iVertex] = latVertex[vertexIndices]
        grid_corner_lon[cellIndices, iVertex] = lonVertex[vertexIndices]

    outFile.variables['grid_corner_lat'][:] = grid_corner_lat[:]
    outFile.variables['grid_corner_lon'][:] = grid_corner_lon[:]

    # Update history attribute of netCDF file
    if hasattr(inFile, 'history'):
        newhist = '\n'.join([getattr(inFile, 'history'),
                             ' '.join(sys.argv[:])])
    else:
        newhist = sys.argv[:]
    setattr(outFile, 'history', newhist)

    inFile.close()
    outFile.close()  # }}}


def lat_lon_file_to_scrip(inFileName, scripFileName, latVarName='lat',
                          lonVarName='lon'):  # {{{
    '''
    Given an MPAS mesh file, create a SCRIP file based on the mesh.

    Parameters
    ----------
    inFileName : str
        The path of the file containing a lat-lon grid

    scripFileName : str
        The path to which the SCRIP file should be written

    latVarName, lonVarName : str, optional
        The name of the latitude and longitude variables in the grid file

    Authors
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/23/2017
    '''
    def interp_extrap_corner(inField):
        outField = numpy.zeros(len(inField)+1)
        outField[1:-1] = 0.5*(inField[0:-1] + inField[1:])
        # extrapolate the ends
        outField[0] = 1.5*inField[0] - 0.5*inField[1]
        outField[-1] = 1.5*inField[-1] - 0.5*inField[-2]
        return outField

    inFile = netCDF4.Dataset(inFileName, 'r')
    outFile = netCDF4.Dataset(scripFileName, 'w')

    # Get info from input file
    lat = numpy.array(inFile.variables[latVarName][:], float)
    lon = numpy.array(inFile.variables[lonVarName][:], float)
    if 'degree' in inFile.variables[latVarName].units:
        units = 'degrees'
    else:
        units = 'radians'

    # interp/extrap corners
    lonCorner = interp_extrap_corner(lon)
    latCorner = interp_extrap_corner(lat)

    _write_lat_lon_scrip(outFile, lat, lon, latCorner, lonCorner, units)

    # Update history attribute of netCDF file
    if hasattr(inFile, 'history'):
        newhist = '\n'.join([getattr(inFile, 'history'),
                             ' '.join(sys.argv[:])])
    else:
        newhist = sys.argv[:]
    setattr(outFile, 'history', newhist)

    inFile.close()
    outFile.close()  # }}}


def lat_lon_array_to_scrip(latCorner, lonCorner, scripFileName,
                           units='degrees'):  # {{{
    '''
    Given an MPAS mesh file, create a SCRIP file based on the mesh.

    Parameters
    ----------
    latCorner,  lonCorner : str
        One dimensional arrays defining the latitude and longitude coordinates
        of grid corners on the destination grid

    scripFileName : str
        The path to which the SCRIP file should be written

    units : {'degrees', 'radians'}, optional
        The units of `latCorner` and `lonCorner`

    Authors
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/23/2017
    '''

    lon = 0.5*(lonCorner[0:-1] + lonCorner[1:])
    lat = 0.5*(latCorner[0:-1] + latCorner[1:])

    outFile = netCDF4.Dataset(scripFileName, 'w')

    _write_lat_lon_scrip(outFile, lat, lon, latCorner, lonCorner, units)

    # Add history attribute to netCDF file
    setattr(outFile, 'history', sys.argv[:])

    outFile.close()  # }}}


def _create_scrip(outFile, grid_size, grid_corners, grid_rank, units):  # {{{
    '''
    Given a SCRIP files, creates common variables and writes common values used
    in various types of SCRIP files.

    Parameters
    ----------
    outFile : file pointer
        A SCRIP file opened in write mode

    grid_size : int
        The number of elements in the grid or mesh

    grid_corners : int
        The number of corners in the grid or mesh

    grid_rank : int
        The dimensionality of the grid (1 for mesh, 2 for grid)

    units : {'degrees', 'radians'}
        The units for latitude and longitude

    Authors
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/20/2017
    '''
    # Write to output file
    # Dimensions
    outFile.createDimension("grid_size", grid_size)
    outFile.createDimension("grid_corners", grid_corners)
    outFile.createDimension("grid_rank", grid_rank)

    # Variables
    grid_center_lat = outFile.createVariable('grid_center_lat', 'f8',
                                             ('grid_size',))
    grid_center_lat.units = units
    grid_center_lon = outFile.createVariable('grid_center_lon', 'f8',
                                             ('grid_size',))
    grid_center_lon.units = units
    grid_corner_lat = outFile.createVariable('grid_corner_lat', 'f8',
                                             ('grid_size', 'grid_corners'))
    grid_corner_lat.units = units
    grid_corner_lon = outFile.createVariable('grid_corner_lon', 'f8',
                                             ('grid_size', 'grid_corners'))
    grid_corner_lon.units = units
    grid_imask = outFile.createVariable('grid_imask', 'i4', ('grid_size',))
    grid_imask.units = 'unitless'
    outFile.createVariable('grid_dims', 'i4', ('grid_rank',))  # }}}


def _write_lat_lon_scrip(outFile, lat, lon, latCorner, lonCorner,
                         units):  # {{{
    '''
    Given a SCRIP files, creates common variables and writes common values used
    in various types of SCRIP files.

    Parameters
    ----------
    outFile : file pointer
        A SCRIP file opened in write mode

    lat, lon : 1D numpy.array
        Latitude and longitude arrays for cell centers on the grid or mesh

    latCorner, lonCorner : 1D numpy.array
        Latitude and longitude arrays for cell corners on the grid or mesh

    units : {'degrees', 'radians'}
        The units for latitude and longitude

    Authors
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    02/23/2017
    '''
    def unwrap_corners(inField):
        outField = numpy.zeros(((inField.shape[0]-1)*(inField.shape[1]-1), 4))
        # corners are counterclockwise
        outField[:, 0] = inField[0:-1, 0:-1].flat
        outField[:, 1] = inField[0:-1, 1:].flat
        outField[:, 2] = inField[1:, 1:].flat
        outField[:, 3] = inField[1:, 0:-1].flat

        return outField

    nLat = len(lat)
    nLon = len(lon)

    grid_size = nLat*nLon

    _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                  grid_rank=2, units=units)

    (Lon, Lat) = numpy.meshgrid(lon, lat)
    (LonCorner, LatCorner) = numpy.meshgrid(lonCorner, latCorner)

    outFile.variables['grid_center_lat'][:] = Lat.flat
    outFile.variables['grid_center_lon'][:] = Lon.flat
    outFile.variables['grid_dims'][:] = [nLon, nLat]
    outFile.variables['grid_imask'][:] = 1

    outFile.variables['grid_corner_lat'][:] = unwrap_corners(LatCorner)
    outFile.variables['grid_corner_lon'][:] = unwrap_corners(LonCorner) # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
