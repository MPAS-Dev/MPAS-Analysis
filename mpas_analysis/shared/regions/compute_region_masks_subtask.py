# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

import os
import xarray as xr
import numpy
import shapely.geometry
import json
from multiprocessing import Pool
import progressbar
from functools import partial
import mpas_tools.conversion

from geometric_features import read_feature_collection, GeometricFeatures
from geometric_features.aggregation import get_aggregator_by_name

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_region_mask
from mpas_analysis.shared.io import write_netcdf


def get_feature_list(geojsonFileName):
    """
    Builds a list of features found in the geojson file
    """
    # Authors
    # -------
    # Xylar Asay-Davis
    featureList = []
    with open(geojsonFileName) as f:
        featureData = json.load(f)

        for feature in featureData['features']:
            name = feature['properties']['name']
            featureList.append(name)
    return featureList


def compute_mpas_region_masks(geojsonFileName, meshFileName, maskFileName,
                              featureList=None, logger=None, processCount=1,
                              chunkSize=1000, showProgress=True,
                              useMpasMaskCreator=False, dir=None):
    """
    Build a region mask file from the given MPAS mesh and geojson file defining
    a set of regions.
    """
    if os.path.exists(maskFileName):
        return

    if useMpasMaskCreator:
        dsMesh = xr.open_dataset(meshFileName)
        fcMask = read_feature_collection(geojsonFileName)
        dsMasks = mpas_tools.conversion.mask(dsMesh=dsMesh, fcMask=fcMask,
                                             logger=logger, dir=dir)

    else:
        with xr.open_dataset(meshFileName) as dsMesh:
            dsMesh = dsMesh[['lonCell', 'latCell']]
            latCell = numpy.rad2deg(dsMesh.latCell.values)

            # transform longitudes to [-180, 180)
            lonCell = numpy.mod(numpy.rad2deg(dsMesh.lonCell.values) + 180.,
                                360.) - 180.

        # create shapely geometry for lonCell and latCell
        cellPoints = [shapely.geometry.Point(x, y) for x, y in
                      zip(lonCell, latCell)]

        regionNames, masks, properties, nChar = compute_region_masks(
            geojsonFileName, cellPoints, maskFileName, featureList, logger,
            processCount, chunkSize, showProgress)

        nCells = len(cellPoints)

        # create a new data array for masks and another for mask names
        if logger is not None:
            logger.info('  Creating and writing masks dataset...')
        nRegions = len(regionNames)
        dsMasks = xr.Dataset()
        dsMasks['regionCellMasks'] = (('nRegions', 'nCells'),
                                      numpy.zeros((nRegions, nCells), dtype=bool))
        dsMasks['regionNames'] = (('nRegions'),
                                  numpy.zeros((nRegions),
                                              dtype='|S{}'.format(nChar)))

        for index in range(nRegions):
            regionName = regionNames[index]
            mask = masks[index]
            dsMasks['regionCellMasks'][index, :] = mask
            dsMasks['regionNames'][index] = regionName

        for propertyName in properties:
            dsMasks[propertyName] = (('nRegions'), properties[propertyName])

    write_netcdf(dsMasks, maskFileName)


def compute_lon_lat_region_masks(geojsonFileName, lon, lat, maskFileName,
                                 featureList=None, logger=None, processCount=1,
                                 chunkSize=1000, showProgress=True,
                                 lonDim='lon', latDim='lat'):
    """
    Build a region mask file from the given lon, lat and geojson file defining
    a set of regions.
    """
    if os.path.exists(maskFileName):
        return

    nLon = len(lon)
    nLat = len(lat)

    # make sure -180 <= lon < 180
    lon = numpy.mod(lon + 180., 360.) - 180.

    if lonDim != latDim:
        lon, lat = numpy.meshgrid(lon, lat)

    # create shapely geometry for lonCell and latCell
    cellPoints = [shapely.geometry.Point(x, y) for x, y in
                  zip(lon.ravel(), lat.ravel())]

    regionNames, masks, properties, nChar = compute_region_masks(
        geojsonFileName, cellPoints, maskFileName, featureList, logger,
        processCount, chunkSize, showProgress)

    # create a new data array for masks and another for mask names
    if logger is not None:
        logger.info('  Creating and writing masks dataset...')
    nRegions = len(regionNames)
    dsMasks = xr.Dataset()
    if lonDim == latDim:
        dsMasks['regionCellMasks'] = (('nRegions', lonDim),
                                      numpy.zeros((nRegions, nLon),
                                                  dtype=bool))
    else:
        dsMasks['regionCellMasks'] = (('nRegions', latDim, lonDim),
                                      numpy.zeros((nRegions, nLat, nLon),
                                                  dtype=bool))

    dsMasks['regionNames'] = (('nRegions'),
                              numpy.zeros((nRegions),
                                          dtype='|S{}'.format(nChar)))

    for index in range(nRegions):
        regionName = regionNames[index]
        mask = masks[index]
        if lonDim == latDim:
            dsMasks['regionCellMasks'][index, :] = mask
        else:
            dsMasks['regionCellMasks'][index, :] = mask.reshape([nLat, nLon])
        dsMasks['regionNames'][index] = regionName

    for propertyName in properties:
        dsMasks['{}Regions'.format(propertyName)] = \
            (('nRegions'), properties[propertyName])

    write_netcdf(dsMasks, maskFileName)


def compute_region_masks(geojsonFileName, cellPoints, maskFileName,
                         featureList=None, logger=None, processCount=1,
                         chunkSize=1000, showProgress=True):
    """
    Build a region mask file from the given mesh and geojson file defining
    a set of regions.
    """
    if os.path.exists(maskFileName):
        return

    if logger is not None:
        logger.info('Creating masks file {}'.format(maskFileName))

    if featureList is None:
        # get a list of features for use by other tasks (e.g. to determine
        # plot names)
        featureList = get_feature_list(geojsonFileName)

    nCells = len(cellPoints)

    with open(geojsonFileName) as f:
        featureData = json.load(f)

    regionNames = []
    for feature in featureData['features']:
        name = feature['properties']['name']
        if name not in featureList:
            continue
        regionNames.append(name)

    propertyNames = set()
    for feature in featureData['features']:
        for propertyName in feature['properties']:
            if propertyName not in ['name', 'author', 'tags', 'component',
                                    'object']:
                propertyNames.add(propertyName)

    properties = {}
    for propertyName in propertyNames:
        properties[propertyName] = []
        for feature in featureData['features']:
            if propertyName in feature['properties']:
                propertyVal = feature['properties'][propertyName]
                properties[propertyName].append(propertyVal)
            else:
                properties[propertyName].append('')

    if logger is not None:
        logger.info('  Computing masks from {}...'.format(geojsonFileName))

    masks = []

    nChar = 0
    for feature in featureData['features']:
        name = feature['properties']['name']
        if name not in featureList:
            continue

        if logger is not None:
            logger.info('      {}'.format(name))

        shape = shapely.geometry.shape(feature['geometry'])
        if processCount == 1:
            mask = _contains(shape, cellPoints)
        else:
            nChunks = int(numpy.ceil(nCells / chunkSize))
            chunks = []
            indices = [0]
            for iChunk in range(nChunks):
                start = iChunk * chunkSize
                end = min((iChunk + 1) * chunkSize, nCells)
                chunks.append(cellPoints[start:end])
                indices.append(end)

            partial_func = partial(_contains, shape)
            pool = Pool(processCount)

            if showProgress:
                widgets = ['  ', progressbar.Percentage(), ' ',
                           progressbar.Bar(), ' ', progressbar.ETA()]
                bar = progressbar.ProgressBar(widgets=widgets,
                                              max_value=nChunks).start()
            else:
                bar = None

            mask = numpy.zeros((nCells,), bool)
            for iChunk, maskChunk in \
                    enumerate(pool.imap(partial_func, chunks)):
                mask[indices[iChunk]:indices[iChunk + 1]] = maskChunk
                if showProgress:
                    bar.update(iChunk + 1)
            if showProgress:
                bar.finish()
            pool.terminate()

        nChar = max(nChar, len(name))

        masks.append(mask)

    return regionNames, masks, properties, nChar  # }}}


def _contains(shape, cellPoints):
    mask = numpy.array([shape.contains(point) for point in cellPoints],
                       dtype=bool)
    return mask


class ComputeRegionMasksSubtask(AnalysisTask):  # {{{
    """
    An analysis tasks for computing cell masks for regions defined by geojson
    features

    Attributes
    ----------
    regionGroup : str
        The name of one of the supported region groups (see
        :py:func:`geometric_features.aggregation.get_region_by_name()`)

    aggregationFunction : callable
        An aggregation function returned by
        :py:func:`geometric_features.aggregation.get_region_by_name()`

    geojsonFileName : str
        A geojson file, typically from the MPAS ``geometric_features``
        repository, defining the shapes to be masked

    outFileSuffix : str
        The suffix for the resulting mask file

    featureList : list of str
        A list of features to include or ``None`` for all features

    maskFileName : str
        The name of the output mask file

    obsFileName : str
        The name of an observations file to create masks for.  But default,
        lon/lat are taken from an MPAS restart file

    lonVar, latVar : str
        The name of the longitude and latitude variables in ``obsFileName``

    meshName : str
        The name of the mesh or grid, used as part of the mask file name.
        Default is the MPAS mesh name
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, meshName, subprocessCount=1,
                 obsFileName=None, lonVar='lon', latVar='lat',
                 useMpasMaskCreator=False):
        # {{{
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            The name of one of the supported region groups (see
            :py:func:`geometric_features.aggregation.get_region_by_name()`)

        meshName : str
            The name of the mesh or grid, used as part of the mask file name.
            Default is the MPAS mesh name


        subprocessCount : int, optional
            The number of processes that can be used to make the mask

        obsFileName : str, optional
            The name of an observations file to create masks for.  But default,
            lon/lat are taken from an MPAS restart file

        lonVar, latVar : str, optional
            The name of the longitude and latitude variables in ``obsFileName``

        useMpasMaskCreator : bool, optional
            If ``True``, the mask creator from ``mpas_tools`` will be used
            to create the mask.  Otherwise, python code is used.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        suffix = regionGroup.replace(' ', '')
        subtaskName = '{}_{}'.format(meshName, suffix)

        # call the constructor from the base class (AnalysisTask)
        super(ComputeRegionMasksSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=[])

        self.regionGroup = regionGroup
        self.featureList = None
        self.subprocessCount = subprocessCount

        self.obsFileName = obsFileName
        self.lonVar = lonVar
        self.latVar = latVar
        self.meshName = meshName
        self.useMpasMaskCreator = useMpasMaskCreator
        self.useMpasMesh = self.obsFileName is None
        self.maskFileName = None

        self.aggregationFunction, prefix, date = get_aggregator_by_name(
            self.regionGroup)
        self.outFileSuffix = '{}{}'.format(prefix, date)
        self.geojsonFileName = \
            get_region_mask(self.config,
                            '{}.geojson'.format(self.outFileSuffix))

        if not self.useMpasMaskCreator:
            # because this uses a Pool, it cannot be launched as a separate
            # process
            self.runDirectly = True

        parentTask.add_subtask(self)

        # }}}

    def make_region_mask(self):
        """
        If the geojson mask file has not already been cached in the diagnostics
        or custom diagnostic directories, it will be created in the analysis
        output's masks directory.
        """
        function = self.aggregationFunction
        filename = self.geojsonFileName
        if not os.path.exists(filename):
            gf = GeometricFeatures()
            fc = function(gf)
            fc.to_geojson(filename)

    def expand_region_names(self, regionNames):
        """
        If ``regionNames`` contains ``'all'``, make sure the geojson file exists
        and then return all the region names found in the file.

        Parameters
        ----------
        regionNames : list
            A list of region names

        Returns
        -------
        regionNames : list
            A list of region names
        """
        if 'all' in regionNames:
            self.make_region_mask()
            regionNames = get_feature_list(self.geojsonFileName)
        return regionNames

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

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

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ComputeRegionMasksSubtask, self).setup_and_check()

        if self.useMpasMesh:
            try:
                self.obsFileName = self.runStreams.readpath('restart')[0]
            except ValueError:
                raise IOError('No MPAS restart file found: need at least one '
                              'restart file to perform region masking.')

        maskSubdirectory = build_config_full_path(self.config, 'output',
                                                  'maskSubdirectory')
        make_directories(maskSubdirectory)

        self.maskFileName = get_region_mask(
            self.config, '{}_{}.nc'.format(self.meshName, self.outFileSuffix))

        if not os.path.exists(self.maskFileName):
            # no cached mask file, so let's see if there's already one in the
            # masks subfolder of the output directory

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            self.maskFileName = '{}/{}_{}.nc'.format(maskSubdirectory,
                                                     self.meshName,
                                                     self.outFileSuffix)

        if os.path.exists(self.maskFileName):
            # nothing to do so don't block a bunch of other processes
            self.subprocessCount = 1
        # }}}

    def run_task(self):  # {{{
        """
        Compute the requested climatologies
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if os.path.exists(self.maskFileName):
            return

        # make the geojson file if it doesn't exist
        self.make_region_mask()

        if self.featureList is None:
            # get a list of features for use by other tasks (e.g. to determine
            # plot names)
            self.featureList = get_feature_list(self.geojsonFileName)

        if self.useMpasMesh:

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            make_directories(maskSubdirectory)

            compute_mpas_region_masks(
                self.geojsonFileName, self.obsFileName, self.maskFileName,
                self.featureList, self.logger, self.subprocessCount,
                showProgress=False, useMpasMaskCreator=self.useMpasMaskCreator,
                dir=maskSubdirectory)
        else:

            dsGrid = xr.open_dataset(self.obsFileName)
            latVar = dsGrid[self.latVar]
            lonVar = dsGrid[self.lonVar]
            if len(latVar.dims) > 1 or len(lonVar.dims) > 1:
                raise ValueError('Masking does not support multidimensional'
                                 'lat/lon with dims {}'.format(latVar.dims))

            latDim = latVar.dims[0]
            lonDim = lonVar.dims[0]
            lat = latVar.values
            lon = lonVar.values

            compute_lon_lat_region_masks(
                self.geojsonFileName, lon, lat, self.maskFileName,
                self.featureList, self.logger, self.subprocessCount,
                showProgress=False, lonDim=lonDim, latDim=latDim)

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
