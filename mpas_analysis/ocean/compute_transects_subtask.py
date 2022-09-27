# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
import numpy
import xarray as xr
import os
from collections import OrderedDict

from pyremap import PointCollectionDescriptor

from mpas_tools.viz import mesh_to_triangles
from mpas_tools.transects import subdivide_great_circle, \
    cartesian_to_great_circle_distance
from mpas_tools.viz.transects import find_transect_cells_and_weights, \
    make_triangle_tree
from mpas_tools.ocean.transects import find_transect_levels_and_weights, \
    interp_mpas_to_transect_triangle_nodes, \
    interp_transect_grid_to_transect_triangle_nodes

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    get_climatology_op_directory

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.ocean.utility import compute_zmid

from mpas_analysis.shared.interpolation import interp_1d


class ComputeTransectsSubtask(RemapMpasClimatologySubtask):
    """
    A subtask for remapping climatologies to transect points

    Attributes
    ----------

    obsDatasets : TransectsObservations
        A dictionary of observational datasets

    verticalComparisonGridName : {'obs', 'mpas'} or any str
        The vertical grid name on which to compare MPAS data with
        observations. 'obs' indicates the locations of the original
        observations; 'mpas' is the vertical locations of MPAS points,
        remapped to the observation latitude/longitude. If any other,
        string, verticalComparisonGrid should be a 1D numpy array and this
        name should be a useful (and unique) description of that grid.

    verticalComparisonGrid : 1D numpy array
        The vertical grid on which to compare MPAS data with observations
        if ``verticalComparisonGridName`` is not 'obs' or 'mpas'.  The
        values should be elevations (in m, typically negative).

    transectNumber : ``xarray.DataArray``
        For each point in the point collection after remapping, the index of
        the transect it belongs to (so that remapped results can be separated
        back into individual transects for plotting)

    transectCollectionName : str
        A name that describes the collection of transects (e.g. the name
        of the collection of observations) used to name the
        destination "mesh" for regridding

    collectionDescriptor : ``PointCollectionDescriptor``
        The mesh descriptor for the collection of all points in all transects,
        used for remapping

    zMid : ``xarray.DataArray``
        Vertical coordinate at the center of layers, used to interpolate to
        reference depths
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 transectCollectionName, variableList, seasons, obsDatasets,
                 verticalComparisonGridName='obs', verticalComparisonGrid=None,
                 subtaskName='remapTransects'):

        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced a climatology to be remapped and plotted
            as a transect

        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding

        variableList : list of str
            A list of variable names in ``timeSeriesStatsMonthly`` to be
            included in the climatologies

        seasons : list of str
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        obsDatasets : TransectsObservations
            A dictionary of observational datasets

        verticalComparisonGridName : {'obs', 'mpas'} or any str, optional
            The vertical grid name on which to compare MPAS data with
            observations. 'obs' indicates the locations of the original
            observations; 'mpas' is the vertical locations of MPAS points,
            remapped to the observation latitude/longitude. If any other,
            string, verticalComparisonGrid should be a 1D numpy array and this
            name should be a useful (and unique) description of that grid.

        verticalComparisonGrid : 1D numpy array, optional
            The vertical grid on which to compare MPAS data with observations
            if ``verticalComparisonGridName`` is not 'obs' or 'mpas'.  The
            values should be elevations (in m, typically negative).

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super(ComputeTransectsSubtask, self).__init__(
            mpasClimatologyTask, parentTask,
            climatologyName=climatologyName, variableList=variableList,
            seasons=seasons, subtaskName=subtaskName)

        self.obsDatasets = obsDatasets
        self.transectCollectionName = transectCollectionName
        self.verticalComparisonGridName = verticalComparisonGridName
        self.verticalComparisonGrid = verticalComparisonGrid
        self.transectNumber = None
        self.x = None
        self.collectionDescriptor = None
        self.maxLevelCell = None
        self.zMid = None
        self.remap = self.obsDatasets.horizontalResolution != 'mpas'
        if self.obsDatasets.horizontalResolution == 'mpas' and \
                self.verticalComparisonGridName != 'mpas':
            raise ValueError('If the horizontal comparison grid is "mpas", the '
                             'vertical grid must also be "mpas".')

    def setup_and_check(self):
        """
        Creates a PointCollectionDescriptor describing all the points in the
        transects to remap to.  Keeps track of which transects index each point
        belongs to.

        Raises
        ------
        IOError :
            If a restart file is not available from which to read mesh
            information or if no history files are available from which to
            compute the climatology in the desired time range.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.remap:
            transectNumber = []
            lats = []
            lons = []
            x = []
            obsDatasets = self.obsDatasets.get_observations()
            datasets = list(obsDatasets.values())
            for transectIndex, ds in enumerate(datasets):
                localLats = list(ds.lat.values)
                localLons = list(ds.lon.values)
                localX = list(ds.x.values)
                localIndices = [transectIndex for _ in localLats]
                lats.extend(localLats)
                lons.extend(localLons)
                x.extend(localX)
                transectNumber.extend(localIndices)

            self.transectNumber = xr.DataArray.from_dict(
                {'dims': ('nPoints',),
                 'data': transectNumber})

            self.x = xr.DataArray.from_dict(
                {'dims': ('nPoints',),
                 'data': x})

            self.collectionDescriptor = PointCollectionDescriptor(
                lats, lons, collectionName=self.transectCollectionName,
                units='degrees', outDimension='nPoints')

            self.add_comparison_grid_descriptor(self.transectCollectionName,
                                                self.collectionDescriptor)

            for transectName in obsDatasets:
                obsDatasets[transectName].close()

        # then, call setup_and_check from the parent class
        # (RemapMpasClimatologySubtask)
        super().setup_and_check()

    def run_task(self):
        """
        Compute climatologies of melt rates from E3SM/MPAS output

        This function has been overridden to compute  ``zMid`` based on data
        from a restart file for later use in vertically interpolating to
        reference depths.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, get maxLevelCell and zMid, needed for masking

        dsMesh = xr.open_dataset(self.restartFileName)
        dsMesh = dsMesh.isel(Time=0)

        self.maxLevelCell = dsMesh.maxLevelCell - 1

        if self.remap:
            zMid = compute_zmid(dsMesh.bottomDepth, dsMesh.maxLevelCell-1,
                                dsMesh.layerThickness)

            self.zMid = \
                xr.DataArray.from_dict({'dims': ('nCells', 'nVertLevels'),
                                        'data': zMid})

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform masking and possibly horizontal remapping
        super(ComputeTransectsSubtask, self).run_task()

        obsDatasets = self.obsDatasets.get_observations()

        if self.remap:

            self.logger.info('Interpolating each transect vertically...')
            # vertically interpolate and write out each transect
            for season in self.seasons:

                remappedFileName = self.get_remapped_file_name(
                    season, comparisonGridName=self.transectCollectionName)

                ds = xr.open_dataset(remappedFileName, decode_times=False)
                transectNames = list(obsDatasets.keys())
                for transectIndex, transectName in enumerate(transectNames):
                    self.logger.info('  {}'.format(transectName))
                    dsObs = obsDatasets[transectName]
                    outFileName = self.get_remapped_file_name(
                        season, comparisonGridName=transectName)
                    outObsFileName = self.obsDatasets.get_out_file_name(
                        transectName, self.verticalComparisonGridName)
                    self._vertical_interp(ds, transectIndex, dsObs,
                                          outFileName, outObsFileName)
                ds.close()

            for transectName in obsDatasets:
                obsDatasets[transectName].close()

        else:
            self._compute_mpas_transects(dsMesh)

        dsMesh.close()

    def customize_masked_climatology(self, climatology, season):
        """
        Add zMid to the climatologies

        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
            the climatology data set

        season : str
            The name of the season to be masked

        Returns
        -------
        climatology : ``xarray.Dataset`` object
            the modified climatology data set
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        zIndex = xr.DataArray.from_dict(
            {'dims': ('nVertLevels',),
             'data': numpy.arange(climatology.sizes['nVertLevels'])})

        cellMask = zIndex <= self.maxLevelCell

        for variableName in self.variableList:
            climatology[variableName] = \
                climatology[variableName].where(cellMask)

        if self.remap:
            climatology['zMid'] = self.zMid

        climatology = climatology.transpose('nVertLevels', 'nCells')

        return climatology

    def customize_remapped_climatology(self, climatology, comparisonGridNames,
                                       season):
        """
        Add the transect index to the data set

        Parameters
        ----------
        climatology : ``xarray.Dataset```
            The MPAS climatology data set that has been remapped

        comparisonGridNames : {'latlon', 'antarctic'}
            The name of the comparison grid to use for remapping.

        season : str
            The name of the season to be masked

        Returns
        -------
        climatology : ``xarray.Dataset```
            The same data set with any custom fields added or modifications
            made
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        climatology['transectNumber'] = self.transectNumber

        climatology['x'] = self.x

        if 'nCells' in climatology.dims:
            climatology = climatology.rename({'nCells': 'nPoints'})

        dims = ['nPoints', 'nVertLevels']
        if 'nv' in climatology.dims:
            dims.append('nv')
        climatology = climatology.transpose(*dims)

        return climatology

    def _vertical_interp(self, ds, transectIndex, dsObs, outFileName,
                         outObsFileName):
        """
        Vertically interpolate a transect and write it to a unique file

        Parameters
        ----------
        ds : ``xarray.Dataset``
            The data set containing all transects before vertical interpolation

        transectIndex : int
            The index of the transect to extract

        dsObs : ``xarray.Dataset``
            The obs dataset used if verticalComparisonGridName is 'obs'

        outFileName : str
            The name of the file to which the resulting data set should be
            written

        outObsFileName : str
            The name of the file to which the resulting obs data set should be
            written if it is interpolated
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if os.path.exists(outFileName):
            return

        ds = ds.where(ds.transectNumber == transectIndex, drop=True)

        if self.verticalComparisonGridName == 'mpas':
            z = ds.zMid
            z = z.rename({'nVertLevels': 'nzOut'})
        elif self.verticalComparisonGridName == 'obs':
            z = dsObs.z
            z = z.rename({'nz': 'nzOut'})
        else:
            # a defined vertical grid
            z = (('nzOut', ), self.verticalComparisonGrid)

        if self.verticalComparisonGridName == 'mpas':
            ds = ds.rename({'zMid': 'z', 'nVertLevels': 'nz'})
        else:
            ds['z'] = z
            # remap each variable
            ds = interp_1d(ds, inInterpDim='nVertLevels', inInterpCoord='zMid',
                           outInterpDim='nzOut', outInterpCoord='z')
            ds = ds.rename({'nzOut': 'nz'})

        if self.verticalComparisonGridName != 'obs' and 'nz' in dsObs.dims:
            dsObs['zOut'] = z
            # remap each variable
            dsObs = interp_1d(dsObs, inInterpDim='nz', inInterpCoord='z',
                              outInterpDim='nzOut', outInterpCoord='zOut')
            dsObs = dsObs.rename({'nzOut': 'nz'})
            write_netcdf(dsObs, outObsFileName)

        ds = ds.drop_vars(['validMask', 'transectNumber'])
        write_netcdf(ds, outFileName)

    def get_mpas_transect_file_name(self, transectName):
        """Get the file name for a masked MPAS transect info"""
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config
        mpasMeshName = config.get('input', 'mpasMeshName')

        climatologyOpDirectory = get_climatology_op_directory(config, 'avg')

        comparisonFullMeshName = transectName.replace(' ', '_')

        stageDirectory = '{}/remapped'.format(climatologyOpDirectory)

        directory = '{}/{}_{}_to_{}'.format(
            stageDirectory, self.climatologyName, mpasMeshName,
            comparisonFullMeshName)

        make_directories(directory)

        fileName = '{}/mpas_transect_info.nc'.format(directory)

        return fileName

    def _compute_mpas_transects(self, dsMesh):

        # see if all transects have already been computed
        allExist = True
        transectNames = list(self.obsDatasets.obsFileNames.keys())
        for transectName in transectNames:
            transectInfoFileName = self.get_mpas_transect_file_name(
                transectName)
            if not os.path.exists(transectInfoFileName):
                allExist = False
                break
            obsFileName = self.obsDatasets.get_out_file_name(
                transectName, self.verticalComparisonGridName)
            if not os.path.exists(obsFileName):
                allExist = False
                break

        if allExist:
            return

        dsTris = mesh_to_triangles(dsMesh)

        triangleTree = make_triangle_tree(dsTris)

        for transectName in transectNames:
            obsFileName = self.obsDatasets.get_out_file_name(
                transectName, self.verticalComparisonGridName)
            transectInfoFileName = self.get_mpas_transect_file_name(
                transectName)
            if not os.path.exists(obsFileName) or  \
                    not os.path.exists(transectInfoFileName):
                dsTransect = self.obsDatasets.build_observational_dataset(
                    self.obsDatasets.obsFileNames[transectName], transectName)

                dsTransect.load()
                # make sure lat and lon are coordinates
                for coord in ['lon', 'lat']:
                    dsTransect.coords[coord] = dsTransect[coord]

                if 'z' in dsTransect:
                    transectZ = dsTransect.z
                else:
                    transectZ = None

                dsMpasTransect = find_transect_cells_and_weights(
                    dsTransect.lon, dsTransect.lat, dsTris, dsMesh,
                    triangleTree, degrees=True)

                dsMpasTransect = find_transect_levels_and_weights(
                    dsMpasTransect, dsMesh.layerThickness,
                    dsMesh.bottomDepth, dsMesh.maxLevelCell - 1,
                    transectZ)

                if 'landIceFraction' in dsMesh:
                    interpCellIndices = dsMpasTransect.interpHorizCellIndices
                    interpCellWeights = dsMpasTransect.interpHorizCellWeights
                    landIceFraction = dsMesh.landIceFraction.isel(
                        nCells=interpCellIndices)
                    landIceFraction = (landIceFraction * interpCellWeights).sum(
                        dim='nHorizWeights')
                    dsMpasTransect['landIceFraction'] = landIceFraction

                # use to_netcdf rather than write_netcdf because integer indices
                # are getting converted to floats when xarray reads them back
                # because of _FillValue
                dsMpasTransect.to_netcdf(transectInfoFileName)

                dsTransectOnMpas = xr.Dataset(dsMpasTransect)
                dsTransectOnMpas['x'] = dsMpasTransect.dNode.isel(
                    nSegments=dsMpasTransect.segmentIndices,
                    nHorizBounds=dsMpasTransect.nodeHorizBoundsIndices)

                dsTransectOnMpas['z'] = dsMpasTransect.zTransectNode

                for var in dsTransect.data_vars:
                    dims = dsTransect[var].dims
                    if 'nPoints' in dims and 'nz' in dims:
                        da = dsTransect[var]
                        da = self._interp_obs_to_mpas(da, dsMpasTransect)
                        dsTransectOnMpas[var] = da

                dsTransectOnMpas.to_netcdf(obsFileName)

        for transectName in transectNames:
            transectInfoFileName = self.get_mpas_transect_file_name(
                transectName)
            dsMpasTransect = xr.open_dataset(transectInfoFileName)

            for season in self.seasons:
                maskedFileName = self.get_masked_file_name(season)
                with xr.open_dataset(maskedFileName) as dsMask:
                    dsOnMpas = xr.Dataset(dsMpasTransect)
                    for var in dsMask.data_vars:
                        dims = dsMask[var].dims
                        if 'nCells' in dims and 'nVertLevels' in dims:
                            dsOnMpas[var] = \
                                interp_mpas_to_transect_triangle_nodes(
                                    dsMpasTransect, dsMask[var])

                    outFileName = self.get_remapped_file_name(
                        season, comparisonGridName=transectName)
                    dsOnMpas.to_netcdf(outFileName)

    def _interp_obs_to_mpas(self, da, dsMpasTransect, threshold=0.1):
        """
        Interpolate observations to the native MPAS transect with masking
        """
        daMask = da.notnull()
        da = da.where(daMask, 0.)
        da = interp_transect_grid_to_transect_triangle_nodes(
            dsMpasTransect, da)
        daMask = interp_transect_grid_to_transect_triangle_nodes(
            dsMpasTransect, daMask)
        da = (da / daMask).where(daMask > threshold)
        return da


class TransectsObservations(object):
    """
    A class for loading and manipulating transect observations

    Attributes
    ----------

    config : mpas_tools.config.MpasConfigParser
        Configuration options

    obsFileNames : OrderedDict
        The names of transects and the file names of the corresponding
        observations for a transect

    horizontalResolution : str
        'obs' for the obs as they are, 'mpas' for the native MPAS mesh, or a
        size in km if subdivision of the observational transect is desired.

    transectCollectionName : str
        A name that describes the collection of transects (e.g. the name
        of the collection of observations) used to name the
        destination "mesh" for regridding

    obsDatasets : OrderedDict
        A dictionary of observational datasets
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, obsFileNames, horizontalResolution,
                 transectCollectionName):
        """
        Construct the object, setting the observations dictionary to None.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        obsFileNames : OrderedDict
            The names of transects and the file names of the corresponding
            observations for a transect

        horizontalResolution : str
            'obs' for the obs as they are, 'mpas' for the native MPAS mesh, or a
            size in km if subdivision of the observational transect is desired.

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.obsDatasets = None
        self.config = config
        self.obsFileNames = obsFileNames
        if horizontalResolution not in ['obs', 'mpas']:
            horizontalResolution = float(horizontalResolution)
        self.horizontalResolution = horizontalResolution
        self.transectCollectionName = transectCollectionName

    def get_observations(self):

        """
        Read in and set up the observations.

        Returns
        -------
        obsDatasets : OrderedDict
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.horizontalResolution == 'mpas':
            # by default, we don't do anything for transects on the native grid
            # but subclasses might need to do something
            return None

        obsDatasets = OrderedDict()
        for name in self.obsFileNames:
            outFileName = self.get_out_file_name(name)
            if os.path.exists(outFileName):
                dsObs = xr.open_dataset(outFileName)
                dsObs.load()
            else:
                dsObs = self.build_observational_dataset(
                    self.obsFileNames[name], name)

                dsObs.load()
                # make sure lat and lon are coordinates
                for coord in ['lon', 'lat']:
                    dsObs.coords[coord] = dsObs[coord]

                dsObs = self._add_distance(dsObs)
                write_netcdf(dsObs, outFileName)
            obsDatasets[name] = dsObs

        return obsDatasets

    def build_observational_dataset(self, fileName, transectName):
        """
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        transectName : str
            transect name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)

        # observations are expected to have horizontal dimension nPoints and
        # vertical dimension nz, as well as horizontal coordinates lat and lon
        # and vertical coordinate z.  Override this function if these need to
        # be renamed from the observations file.

        return dsObs

    def get_out_file_name(self, transectName,
                          verticalComparisonGridName='obs'):
        """
        Given config options, the name of a field and a string identifying the
        months in a seasonal climatology, returns the full path for MPAS
        climatology files before and after remapping.

        Parameters
        ----------
        transectName : str
            The name of the transect

        verticalComparisonGridName : {'obs', 'mpas'} or any str, optional
            The vertical grid name on which to compare MPAS data with
            observations. 'obs' indicates the locations of the original
            observations; 'mpas' is the vertical locations of MPAS points,
            remapped to the observation latitude/longitude. If any other,
            string, verticalComparisonGrid should be a 1D numpy array and this
            name should be a useful (and unique) description of that grid.

        Returns
        -------
        fileName : str
            The path to the climatology file for the specified season.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        remappedDirectory = build_config_full_path(
            config=config, section='output',
            relativePathOption='remappedClimSubdirectory',
            relativePathSection='oceanObservations')

        make_directories(remappedDirectory)

        transectSuffix = transectName.replace(' ', '_')

        if verticalComparisonGridName == 'obs':
            fileName = '{}/{}_{}.nc'.format(
                remappedDirectory, self.transectCollectionName, transectSuffix)
        else:
            fileName = '{}/{}_{}_{}.nc'.format(
                remappedDirectory, self.transectCollectionName, transectSuffix,
                verticalComparisonGridName)
        return fileName

    def _add_distance(self, dsObs):
        """
        Add a distance coordinate for the transect.  If a horizontal resolution
        for subdivision is provided, subdivide each segment of the transect so
        the horizontal resolution is at least as high as the requested
        resolution
        """

        lat = numpy.deg2rad(dsObs.lat.values)
        lon = numpy.deg2rad(dsObs.lon.values)

        earth_radius = 6.371e6  # Radius of earth in meters

        x = earth_radius * numpy.cos(lat) * numpy.cos(lon)
        y = earth_radius * numpy.cos(lat) * numpy.sin(lon)
        z = earth_radius * numpy.sin(lat)

        if self.horizontalResolution == 'obs':
            dIn = cartesian_to_great_circle_distance(x, y, z, earth_radius)
            dsObs['x'] = (('nPoints',), dIn)
        elif self.horizontalResolution != 'mpas':
            # subdivide
            xOut, yOut, zOut, dIn, dOut = subdivide_great_circle(
                x, y, z, 1e3*self.horizontalResolution, earth_radius)

            dsObs['xIn'] = (('nPoints',), dIn)
            dsObs['xOut'] = (('nPointsOut',), dOut)

            # interpolate fields without and with vertical dimension
            dsObs = interp_1d(dsObs, inInterpDim='nPoints',
                              inInterpCoord='xIn', outInterpDim='nPointsOut',
                              outInterpCoord='xOut')
            dsObs = dsObs.drop_vars(['xIn'])
            dsObs = dsObs.rename({'nPointsOut': 'nPoints', 'xOut': 'x'})

        return dsObs
