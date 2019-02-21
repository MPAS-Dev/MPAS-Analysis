# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy
import xarray as xr
import os
from collections import OrderedDict

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.shared.mpas_xarray import mpas_xarray

from mpas_analysis.shared.grid import PointCollectionDescriptor

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.ocean.utility import compute_zmid

from mpas_analysis.shared.interpolation import interp_1d


class ComputeTransectsSubtask(RemapMpasClimatologySubtask):  # {{{
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

        # {{{
        '''
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
        '''
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

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Creates a PointCollectionDescriptor describing all the points in the
        transects to remap to.  Keeps track of which transects index each point
        belongs to.

        Raises
        ------
        IOError :
            If a restart file is not available from which to read mesh
            information or if no history files are available from which to
            compute the climatology in the desired time range.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

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
            localIndices = [transectIndex for lat in localLats]
            lats.extend(localLats)
            lons.extend(localLons)
            x.extend(localX)
            transectNumber.extend(localIndices)

        self.transectNumber = xr.DataArray.from_dict(
            {'dims': ('nPoints'),
             'data': transectNumber})

        self.x = xr.DataArray.from_dict(
            {'dims': ('nPoints'),
             'data': x})

        self.collectionDescriptor = PointCollectionDescriptor(
            lats, lons, collectionName=self.transectCollectionName,
            units='degrees', outDimension='nPoints')

        self.add_comparison_grid_descriptor(self.transectCollectionName,
                                            self.collectionDescriptor)

        # then, call setup_and_check from the base class
        # (RemapMpasClimatologySubtask)
        super(ComputeTransectsSubtask, self).setup_and_check()

        for transectName in obsDatasets:
            obsDatasets[transectName].close()

    def run_task(self):  # {{{
        '''
        Compute climatologies of melt rates from E3SM/MPAS output

        This function has been overridden to compute  ``zMid`` based on data
        from a restart file for later use in vertically interpolating to
        reference depths.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, compute zMid and cell mask from the restart file
        with xr.open_dataset(self.restartFileName) as ds:
            ds = mpas_xarray.subset_variables(ds, ['maxLevelCell',
                                                   'bottomDepth',
                                                   'layerThickness'])
            ds = ds.isel(Time=0)

            self.maxLevelCell = ds.maxLevelCell - 1

            zMid = compute_zmid(ds.bottomDepth, ds.maxLevelCell,
                                ds.layerThickness)

            self.zMid = \
                xr.DataArray.from_dict({'dims': ('nCells', 'nVertLevels'),
                                        'data': zMid})
            ds.close()

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform the horizontal remapping
        super(ComputeTransectsSubtask, self).run_task()

        obsDatasets = self.obsDatasets.get_observations()

        self.logger.info('Interpolating each transect vertically...')
        # finally, vertically interpolate and write out each transect
        for season in self.seasons:

            remappedFileName = self.get_remapped_file_name(
                season, comparisonGridName=self.transectCollectionName)

            with xr.open_dataset(remappedFileName) as ds:
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

        # }}}

    def customize_masked_climatology(self, climatology, season):  # {{{
        '''
        Add zMid to the climatologys

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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        zIndex = xr.DataArray.from_dict(
            {'dims': ('nVertLevels',),
             'data': numpy.arange(climatology.sizes['nVertLevels'])})

        cellMask = zIndex < self.maxLevelCell

        for variableName in self.variableList:
            climatology[variableName] = \
                climatology[variableName].where(cellMask)

        climatology['zMid'] = self.zMid

        return climatology  # }}}

    def customize_remapped_climatology(self, climatology, comparisonGridNames,
                                       season):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        climatology['transectNumber'] = self.transectNumber

        climatology['x'] = self.x

        return climatology  # }}}

    def _vertical_interp(self, ds, transectIndex, dsObs, outFileName,
                         outObsFileName):
        '''
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
        '''
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

        ds = ds.drop(['validMask', 'transectNumber'])
        write_netcdf(ds, outFileName)  # }}}

    # }}}


class TransectsObservations(object):  # {{{
    """
    A class for loading and manipulating transect observations

    Attributes
    ----------

    config :  ``MpasAnalysisConfigParser``
        Configuration options

    obsFileNames : OrderedDict
        The names of transects and the file names of the corresponding
        observations for a transect

    horizontalResolution : str
        'obs' for the obs as they are or a size in km if subdivision is
        desired.

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
                 transectCollectionName):  # {{{
        '''
        Construct the object, setting the observations dictionary to None.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        obsFileNames : OrderedDict
            The names of transects and the file names of the corresponding
            observations for a transect

        horizontalResolution : str
            'obs' for the obs as they are or a size in km if subdivision is
            desired.

        transectCollectionName : str
            A name that describes the collection of transects (e.g. the name
            of the collection of observations) used to name the
            destination "mesh" for regridding
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        self.obsDatasets = None
        self.config = config
        self.obsFileNames = obsFileNames
        if horizontalResolution != 'obs':
            horizontalResolution = float(horizontalResolution)
        self.horizontalResolution = horizontalResolution
        self.transectCollectionName = transectCollectionName

    def get_observations(self):
        # {{{
        '''
        Read in and set up the observations.

        Returns
        -------
        obsDatasets : OrderedDict
            The observational dataset
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

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

                if self.horizontalResolution == 'obs':
                    dsObs = self._add_distance(dsObs)
                else:
                    dsObs = self._subdivide_observations(dsObs)
                write_netcdf(dsObs, outFileName)
            obsDatasets[name] = dsObs

        return obsDatasets  # }}}

    def build_observational_dataset(self, fileName, transectName):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)

        # observations are expected to have horizontal dimension nPoints and
        # vertical dimension nz, as well as horizontal coordinates lat and lon
        # and vertical coordinate z.  Override this function if these need to
        # be renamed from the observations file.

        return dsObs  # }}}

    def get_out_file_name(self, transectName,
                          verticalComparisonGridName='obs'):  # {{{
        '''
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        remappedDirectory = build_config_full_path(
            config=config, section='output',
            relativePathOption='remappedClimSubdirectory',
            relativePathSection='oceanObservations')

        make_directories(remappedDirectory)

        if verticalComparisonGridName == 'obs':
            fileName = '{}/{}_{}.nc'.format(
                remappedDirectory, self.transectCollectionName, transectName)
        else:
            fileName = '{}/{}_{}_{}.nc'.format(
                remappedDirectory, self.transectCollectionName, transectName,
                verticalComparisonGridName)
        return fileName  # }}}

    def _add_distance(self, dsObs):  # {{{
        '''
        Subdivide each segment of the transect so the horizontal resolution
        approximately matches the requested resolution
        '''

        lat = dsObs.lat.values
        lon = dsObs.lon.values

        # compute the great circle distance between these points
        dxIn = self._haversine(lon[0:-1], lat[0:-1], lon[1:], lat[1:])

        xIn = numpy.zeros(lat.shape)
        xIn[1:] = numpy.cumsum(dxIn)

        dsObs['x'] = (('nPoints',), xIn)
        return dsObs  # }}}

    def _subdivide_observations(self, dsObs):  # {{{
        '''
        Subdivide each segment of the transect so the horizontal resolution
        approximately matches the requested resolution
        '''

        lat = dsObs.lat.values
        lon = dsObs.lon.values

        # compute the great circle distance between these points
        dxIn = self._haversine(lon[0:-1], lat[0:-1], lon[1:], lat[1:])

        nSegments = numpy.maximum(
            (dxIn / self.horizontalResolution + 0.5).astype(int), 1)

        xIn = numpy.zeros(lat.shape)
        xIn[1:] = numpy.cumsum(dxIn)

        outIndex = []
        for index in range(len(xIn) - 1):
            n = nSegments[index]
            outIndex.extend(index + numpy.arange(0, n) / n)
        outIndex.append(len(xIn) - 1)

        xOut = numpy.interp(outIndex, numpy.arange(len(xIn)), xIn)

        dsObs['xIn'] = (('nPoints',), xIn)
        dsObs['xOut'] = (('nPointsOut',), xOut)

        # interpolate fields without and with vertical dimension
        dsObs = interp_1d(dsObs, inInterpDim='nPoints',
                          inInterpCoord='xIn', outInterpDim='nPointsOut',
                          outInterpCoord='xOut')
        dsObs = dsObs.drop(['xIn'])
        dsObs = dsObs.rename({'nPointsOut': 'nPoints', 'xOut': 'x'})
        return dsObs  # }}}

    def _haversine(self, lon1, lat1, lon2, lat2):  # {{{
        """
        Calculate the great circle distance in km between two points on the
        earth (specified in decimal degrees). Based on
        https://stackoverflow.com/a/4913653
        """
        # convert decimal degrees to radians
        lon1 = numpy.deg2rad(lon1)
        lat1 = numpy.deg2rad(lat1)
        lon2 = numpy.deg2rad(lon2)
        lat2 = numpy.deg2rad(lat2)

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = numpy.sin(dlat / 2.)**2 + numpy.cos(lat1) * numpy.cos(lat2) * \
            numpy.sin(dlon / 2.)**2
        c = 2 * numpy.arcsin(numpy.sqrt(a))
        r = 6371  # Radius of earth in kilometers. Use 3956 for miles
        return c * r  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
