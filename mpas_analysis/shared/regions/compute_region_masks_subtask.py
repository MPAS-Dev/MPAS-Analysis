# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE

import os
import xarray as xr
import json

from geometric_features import read_feature_collection, GeometricFeatures
from geometric_features.aggregation import get_aggregator_by_name
import mpas_tools.conversion
from mpas_tools.logging import check_call

from mpas_analysis.shared.analysis_task import AnalysisTask
from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_region_mask
from mpas_analysis.shared.io import write_netcdf_with_fill


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
                              logger=None, processCount=1,
                              multiprocessingMethod='spawn', chunkSize=1000,
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
        write_netcdf_with_fill(dsMasks, maskFileName)

    else:
        args = ['compute_mpas_region_masks',
                '-m', meshFileName,
                '-g', geojsonFileName,
                '-o', maskFileName,
                '-t', 'cell',
                '--chunk_size', '{}'.format(chunkSize),
                '--process_count', '{}'.format(processCount),
                '--multiprocessing_method', '{}'.format(multiprocessingMethod)]
        check_call(args=args, logger=logger)


def compute_lon_lat_region_masks(gridFileName, lonVar, latVar, geojsonFileName,
                                 maskFileName, logger=None, processCount=1,
                                 multiprocessingMethod='spawn', chunkSize=1000):
    """
    Build a region mask file from the given lon, lat and geojson file defining
    a set of regions.
    """
    if os.path.exists(maskFileName):
        return

    args = ['compute_lon_lat_region_masks',
            '-i', gridFileName,
            '--lon', lonVar,
            '--lat', latVar,
            '-g', geojsonFileName,
            '-o', maskFileName,
            '--chunk_size', '{}'.format(chunkSize),
            '--process_count', '{}'.format(processCount),
            '--multiprocessing_method', '{}'.format(multiprocessingMethod)]
    check_call(args=args, logger=logger)


class ComputeRegionMasksSubtask(AnalysisTask):
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

    maskFileName : str
        The name of the output mask file

    obsFileName : str
        The name of an observations file to create masks for.  By default,
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
        self.date = date
        self.outFileSuffix = '{}{}'.format(prefix, date)
        self.geojsonFileName = \
            get_region_mask(self.config,
                            '{}.geojson'.format(self.outFileSuffix))

        if not self.useMpasMaskCreator:
            # because this uses a Pool, it cannot be launched as a separate
            # process
            self.runDirectly = True

        parentTask.add_subtask(self)

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

    def setup_and_check(self):
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
                self.obsFileName = self.runStreams.readpath('mesh')[0]
            except ValueError:
                raise IOError(
                    'The MPAS mesh file could not be found: needed to perform '
                    'region masking.'
                )

        maskSubdirectory = build_config_full_path(self.config, 'output',
                                                  'maskSubdirectory')
        make_directories(maskSubdirectory)

        self.maskFileName = get_region_mask(
            self.config, '{}_{}.nc'.format(self.meshName, self.outFileSuffix))

        if not os.path.exists(self.maskFileName):
            # no cached mask file, so let's see if there's already one in the
            # masks subdirectory of the output directory

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            self.maskFileName = '{}/{}_{}.nc'.format(maskSubdirectory,
                                                     self.meshName,
                                                     self.outFileSuffix)

        if os.path.exists(self.maskFileName):
            # nothing to do so don't block a bunch of other processes
            self.subprocessCount = 1

    def run_task(self):
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

        multiprocessingMethod = self.config.get('execute',
                                                'multiprocessingMethod')

        if self.useMpasMesh:

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            make_directories(maskSubdirectory)

            compute_mpas_region_masks(
                self.geojsonFileName, self.obsFileName, self.maskFileName,
                self.logger, self.subprocessCount,
                multiprocessingMethod=multiprocessingMethod,
                useMpasMaskCreator=self.useMpasMaskCreator,
                dir=maskSubdirectory)
        else:
            compute_lon_lat_region_masks(
                self.obsFileName, self.lonVar, self.latVar,
                self.geojsonFileName, self.maskFileName, self.logger,
                self.subprocessCount,
                multiprocessingMethod=multiprocessingMethod)
