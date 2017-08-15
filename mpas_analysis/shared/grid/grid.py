'''
Classes for describing meshes and grids , including creating SCRIP files,
used to create mapping files

Classes
-------
MpasMeshDescriptor - describes an MPAS mesh

LatLonGridDescriptor - describes a lat-lon grid

ProjectionGridDescriptor - describes a logically rectangular grid on a pyproj
    projection

Author
------
Xylar Asay-Davis
'''

import netCDF4
import numpy
import sys
import pyproj
import xarray


class MeshDescriptor(object):  # {{{
    '''
    A class for describing a mesh

    Author
    ------
    Xylar Asay-Davis
    '''

    def __init__(self):  # {{{
        '''
        Constructor creates a common ``meshName`` member variable, ``None`` by
        default.  Each Subclass should define or use input arguments to set
        ``meshName`` to a short description of the mesh or grid.

        Author
        ------
        Xylar Asay-Davis
        '''

        self.meshName = None  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Subclasses should overload this method to write a SCRIP file based on
        the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written

        Authors
        ------
        Xylar Asay-Davis
        '''

        return  # }}}

    # }}}


class MpasMeshDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing an MPAS mesh

    Author
    ------
    Xylar Asay-Davis
    '''

    def __init__(self, fileName, meshName=None):  # {{{
        '''
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh

        meshName : str, optional
            The name of the MPAS mesh (e.g. ``'oEC60to30'`` or
            ``'oRRS18to6'``).  If not provided, the data set in ``fileName``
            must have a global attribute ``meshName`` that will be used
            instead.

        Author
        ------
        Xylar Asay-Davis
        '''

        ds = xarray.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            self.meshName = ds.attrs['meshName']
        else:
            self.meshName = meshName

        self.fileName = fileName
        self.regional = True

        # build coords
        self.coords = {'latCell': {'dims': 'nCells',
                                   'data': ds.latCell.values,
                                   'attrs': {'units': 'radians'}},
                       'lonCell': {'dims': 'nCells',
                                   'data': ds.lonCell.values,
                                   'attrs': {'units': 'radians'}}}
        self.dims = ['nCells']
        self.dimSize = [ds.dims[dim] for dim in self.dims]
        ds.close()  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Given an MPAS mesh file, create a SCRIP file based on the mesh.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written

        Authors
        ------
        Xylar Asay-Davis
        '''
        self.scripFileName = scripFileName

        inFile = netCDF4.Dataset(self.fileName, 'r')
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
                      grid_rank=1, units='radians', meshName=self.meshName)

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
# }}}


class LatLonGridDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing a lat-lon grid

    Author
    ------
    Xylar Asay-Davis
    '''
    def __init__(self):  # {{{
        '''
        Constructor stores the file name

        Parameters
        ----------
        fileName : str
            The path of the file containing the MPAS mesh

        Author
        ------
        Xylar Asay-Davis
        '''
        self.regional = False
        self.meshName = None  # }}}

    def read(self, fileName, latVarName='lat', lonVarName='lon'):  # {{{
        '''
        Read the lat-lon grid from a file with the given lat/lon var names.

        Parameters
        ----------
        fileName : str
            The path of the file containing the lat-lon grid

        latVarName, lonVarName : str, optional
            The name of the latitude and longitude variables in the grid file

        Author
        ------
        Xylar Asay-Davis
        '''
        ds = xarray.open_dataset(fileName)

        if self.meshName is None and 'meshName' in ds.attrs:
            self.meshName = ds.attrs['meshName']

        # Get info from input file
        self.lat = numpy.array(ds[latVarName].values, float)
        self.lon = numpy.array(ds[lonVarName].values, float)
        if 'degree' in ds[latVarName].units:
            self.units = 'degrees'
        else:
            self.units = 'radians'

        self._set_coords(latVarName, lonVarName, ds[latVarName].dims[0],
                         ds[lonVarName].dims[0])

        # interp/extrap corners
        self.lonCorner = _interp_extrap_corner(self.lon)
        self.latCorner = _interp_extrap_corner(self.lat)

        if 'history' in ds.attrs:
            self.history = '\n'.join([ds.attrs['history'],
                                     ' '.join(sys.argv[:])])
        else:
            self.history = sys.argv[:]  # }}}

    def create(self, latCorner, lonCorner, units='degrees'):  # {{{
        '''
        Create the lat-lon grid with the given arrays and units.

        Parameters
        ----------
        latCorner, lonCorner : 1D numpy.arrays
            One dimensional arrays defining the latitude and longitude
            coordinates of grid corners.

        units : {'degrees', 'radians'}, optional
            The units of `latCorner` and `lonCorner`

        Author
        ------
        Xylar Asay-Davis
        '''

        self.latCorner = latCorner
        self.lonCorner = lonCorner
        self.lon = 0.5*(lonCorner[0:-1] + lonCorner[1:])
        self.lat = 0.5*(latCorner[0:-1] + latCorner[1:])
        self.units = units
        self.history = sys.argv[:]
        self._set_coords('lat', 'lon', 'lat', 'lon')  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Given a lat-lon grid file, create a SCRIP file based on the grid.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written

        Authors
        ------
        Xylar Asay-Davis
        '''
        self.scripFileName = scripFileName

        outFile = netCDF4.Dataset(scripFileName, 'w')

        nLat = len(self.lat)
        nLon = len(self.lon)

        grid_size = nLat*nLon

        _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                      grid_rank=2, units=self.units, meshName=self.meshName)

        (Lon, Lat) = numpy.meshgrid(self.lon, self.lat)
        (LonCorner, LatCorner) = numpy.meshgrid(self.lonCorner, self.latCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nLon, nLat]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = _unwrap_corners(LatCorner)
        outFile.variables['grid_corner_lon'][:] = _unwrap_corners(LonCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()  # }}}

    def _set_coords(self, latVarName, lonVarName, latDimName,
                    lonDimName):  # {{{
        '''
        Set up a coords dict with lat and lon
        '''
        self.latVarName = latVarName
        self.lonVarName = lonVarName
        self.coords = {latVarName: {'dims': latDimName,
                                    'data': self.lat,
                                    'attrs': {'units': self.units}},
                       lonVarName: {'dims': lonDimName,
                                    'data': self.lon,
                                    'attrs': {'units': self.units}}}

        self.dims = [latDimName, lonDimName]
        self.dimSize = [len(self.lat), len(self.lon)]

        # set the name of the grid
        dLat = self.lat[1]-self.lat[0]
        dLon = self.lon[1]-self.lon[0]
        if 'degree' in self.units:
            units = 'degree'
        elif 'rad' in self.units:
            units = 'radian'
        else:
            raise ValueError('Could not figure out units {}'.format(
                self.units))
        if self.meshName is None:
            self.meshName = '{}x{}{}'.format(abs(dLat), abs(dLon), units)
        # }}}


class ProjectionGridDescriptor(MeshDescriptor):  # {{{
    '''
    A class for describing a general logically rectangular grid that can be
    defined by a `pyproj` projection.

    Author
    ------
    Xylar Asay-Davis
    '''

    def __init__(self, projection):  # {{{
        '''
        Constructor stores the projection

        Parameters
        ----------
        projection : pyproj.Proj object
            The projection used to map from grid x-y space to latitude and
            longitude

        Author
        ------
        Xylar Asay-Davis
        '''
        self.projection = projection
        self.latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')
        self.regional = True

    def read(self, fileName, meshName=None, xVarName='x', yVarName='y'):  # {{{
        '''
        Given a grid file with x and y coordinates defining the axes of the
        logically rectangular grid, read in the x and y coordinates and
        interpolate/extrapolate to locate corners.

        Parameters
        ----------
        fileName : str
            The path of the file containing the grid data

        meshName : str, optional
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``).  If not
            provided, the data set in ``fileName`` must have a global
            attribute ``meshName`` that will be used instead.

        xVarName, yVarName : str, optional
            The name of the x and y (in meters) variables in the grid file

        Authors
        ------
        Xylar Asay-Davis
        '''

        ds = xarray.open_dataset(fileName)

        if meshName is None:
            if 'meshName' not in ds.attrs:
                raise ValueError('No meshName provided or found in file.')
            self.meshName = ds.attrs['meshName']
        else:
            self.meshName = meshName

        # Get info from input file
        self.x = numpy.array(ds[xVarName].values, float)
        self.y = numpy.array(ds[yVarName].values, float)

        self._set_coords(xVarName, yVarName, ds[xVarName].dims[0],
                         ds[yVarName].dims[0])

        # interp/extrap corners
        self.xCorner = _interp_extrap_corner(self.x)
        self.yCorner = _interp_extrap_corner(self.y)

        # Update history attribute of netCDF file
        if 'history' in ds.attrs:
            self.history = '\n'.join([ds.attrs['history'],
                                     ' '.join(sys.argv[:])])
        else:
            self.history = sys.argv[:]  # }}}

    def create(self, x, y, meshName):  # {{{
        '''
        Given x and y coordinates defining the axes of the logically
        rectangular grid, save the coordinates interpolate/extrapolate to
        locate corners.

        Parameters
        ----------
        x, y : 1D numpy.arrays
            One dimensional arrays defining the x and y coordinates of grid
            cell centers.

        meshName : str
            The name of the grid (e.g. ``'10km_Antarctic_stereo'``)

        Author
        ------
        Xylar Asay-Davis
        '''

        self.meshName = meshName

        self.x = x
        self.y = y

        self._set_coords('x', 'y', 'x', 'y')

        # interp/extrap corners
        self.xCorner = _interp_extrap_corner(self.x)
        self.yCorner = _interp_extrap_corner(self.y)
        self.history = sys.argv[:]  # }}}

    def to_scrip(self, scripFileName):  # {{{
        '''
        Create a SCRIP file based on the grid and projection.

        Parameters
        ----------
        scripFileName : str
            The path to which the SCRIP file should be written

        Authors
        ------
        Xylar Asay-Davis
        '''
        self.scripFileName = scripFileName

        outFile = netCDF4.Dataset(scripFileName, 'w')

        nx = len(self.x)
        ny = len(self.y)

        grid_size = nx*ny

        _create_scrip(outFile, grid_size=grid_size, grid_corners=4,
                      grid_rank=2, units='degrees', meshName=self.meshName)

        (X, Y) = numpy.meshgrid(self.x, self.y)
        (XCorner, YCorner) = numpy.meshgrid(self.xCorner, self.yCorner)

        (Lat, Lon) = self.project_to_lat_lon(X, Y)
        (LatCorner, LonCorner) = self.project_to_lat_lon(XCorner, YCorner)

        outFile.variables['grid_center_lat'][:] = Lat.flat
        outFile.variables['grid_center_lon'][:] = Lon.flat
        outFile.variables['grid_dims'][:] = [nx, ny]
        outFile.variables['grid_imask'][:] = 1

        outFile.variables['grid_corner_lat'][:] = _unwrap_corners(LatCorner)
        outFile.variables['grid_corner_lon'][:] = _unwrap_corners(LonCorner)

        setattr(outFile, 'history', self.history)

        outFile.close()  # }}}

    def project_to_lat_lon(self, X, Y):  # {{{
        '''
        Given X and Y locations of points in a projection, returns the
        corresponding latitude and longitude of each point.

        Parameters
        ----------
        outFile : file pointer
            A SCRIP file opened in write mode

        X, Y : 1D or 2D numpy.array
            X and y arrays of points in in the projection

        Returns
        -------
        Lat, Lon : numpy.array with same shape as X and Y
            the latitude and longitude in degrees of the points

        Authors
        ------
        Xylar Asay-Davis
        '''

        Lon, Lat = pyproj.transform(self.projection, self.latLonProjection,
                                    X, Y, radians=False)

        return (Lat, Lon)  # }}}

    def _set_coords(self, xVarName, yVarName, xDimName, yDimName):  # {{{
        '''
        Set up a coords dict with x, y, lat and lon
        '''
        (X, Y) = numpy.meshgrid(self.x, self.y)
        (Lat, Lon) = self.project_to_lat_lon(X, Y)

        self.coords = {xVarName: {'dims': xDimName,
                                  'data': self.x,
                                  'attrs': {'units': 'meters'}},
                       yVarName: {'dims': yDimName,
                                  'data': self.y,
                                  'attrs': {'units': 'meters'}},
                       'lat': {'dims': (xDimName, yDimName),
                               'data': Lat,
                               'attrs': {'units': 'degrees'}},
                       'lon': {'dims': (xDimName, yDimName),
                               'data': Lon,
                               'attrs': {'units': 'degrees'}}}

        self.dims = [xDimName, yDimName]
        self.dimSize = [len(self.x), len(self.y)]
        # }}}


def _create_scrip(outFile, grid_size, grid_corners, grid_rank, units,
                  meshName):  # {{{
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

    meshName : str
        The name of the mesh

    Authors
    ------
    Xylar Asay-Davis
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
    outFile.createVariable('grid_dims', 'i4', ('grid_rank',))

    outFile.meshName = meshName
    # }}}


def _interp_extrap_corner(inField):  # {{{
    '''Interpolate/extrapolate a 1D field from grid centers to grid corners'''

    outField = numpy.zeros(len(inField)+1)
    outField[1:-1] = 0.5*(inField[0:-1] + inField[1:])
    # extrapolate the ends
    outField[0] = 1.5*inField[0] - 0.5*inField[1]
    outField[-1] = 1.5*inField[-1] - 0.5*inField[-2]
    return outField  # }}}


def _unwrap_corners(inField):
    '''Turn a 2D array of corners into an array of rectangular mesh elements'''
    outField = numpy.zeros(((inField.shape[0]-1)*(inField.shape[1]-1), 4))
    # corners are counterclockwise
    outField[:, 0] = inField[0:-1, 0:-1].flat
    outField[:, 1] = inField[0:-1, 1:].flat
    outField[:, 2] = inField[1:, 1:].flat
    outField[:, 3] = inField[1:, 0:-1].flat

    return outField

# vim: ai ts=4 sts=4 et sw=4 ft=python
