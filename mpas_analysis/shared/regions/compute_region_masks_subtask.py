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

import os
import xarray as xr
import numpy
import shapely.geometry
import json
from multiprocessing import Pool
import progressbar
from functools import partial

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.mpas_xarray import mpas_xarray


def get_feature_list(geojsonFileName):
    '''
    Builds a list of features found in the geojson file
    '''
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


def compute_region_masks(geojsonFileName, meshFileName, maskFileName,
                         featureList=None, logger=None, processCount=1,
                         chunkSize=1000):
    '''
    Build a region mask file from the given mesh and geojson file defining
    a set of regions.
    '''
    if os.path.exists(maskFileName):
        return

    if logger is not None:
        logger.info('Creating masks file {}'.format(maskFileName))

    if featureList is None:
        # get a list of features for use by other tasks (e.g. to determine
        # plot names)
        featureList = get_feature_list(geojsonFileName)

    with xr.open_dataset(meshFileName) as dsMesh:
        dsMesh = mpas_xarray.subset_variables(dsMesh, ['lonCell', 'latCell'])
        latCell = numpy.rad2deg(dsMesh.latCell.values)

        # transform longitudes to [-180, 180)
        lonCell = numpy.mod(numpy.rad2deg(dsMesh.lonCell.values) + 180.,
                            360.) - 180.

    # create shapely geometry for lonCell and latCell
    cellPoints = [shapely.geometry.Point(x, y) for x, y in
                  zip(lonCell, latCell)]

    nCells = len(cellPoints)

    masks = []
    regionNames = []
    nChar = 0
    if logger is not None:
        logger.info('  Computing masks from {}...'.format(geojsonFileName))
    with open(geojsonFileName) as f:
        featureData = json.load(f)

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

            widgets = ['  ', progressbar.Percentage(), ' ',
                       progressbar.Bar(), ' ', progressbar.ETA()]
            bar = progressbar.ProgressBar(widgets=widgets,
                                          maxval=nChunks).start()

            mask = numpy.zeros((nCells,), bool)
            for iChunk, maskChunk in \
                    enumerate(pool.imap(partial_func, chunks)):
                mask[indices[iChunk]:indices[iChunk + 1]] = maskChunk
                bar.update(iChunk + 1)
            bar.finish()
            pool.terminate()

        nChar = max(nChar, len(name))

        masks.append(mask)
        regionNames.append(name)

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

    write_netcdf(dsMasks, maskFileName)

    # }}}


def _contains(shape, cellPoints):
    mask = numpy.array([shape.contains(point) for point in cellPoints],
                       dtype=bool)
    return mask


class ComputeRegionMasksSubtask(AnalysisTask):  # {{{
    '''
    An analysis tasks for computing climatologies from output from the
    ``timeSeriesStatsMonthly`` analysis member.

    Attributes
    ----------
    geojsonFileName : str
        A geojson file, typically from the MPAS ``geometric_features``
        repository, defining the shapes to be masked

    outFileSuffix : str
        The suffix for the resulting mask file

    featureList : list of str
        A list of features to include or ``None`` for all features

    maskFileName : str
        The name of the output mask file

    maskExists : bool
        Whether the mask file already exists

    '''
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, geojsonFileName, outFileSuffix,
                 featureList=None, subtaskName='computeRegionMasks'):
        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        geojsonFileName : str
            A geojson file, typically from the MPAS ``geometric_features``
            repository, defining the shapes to be masked

        outFileSuffix : str
            The suffix for the resulting mask file

        featureList : list of str, optional
            A list of features to include.  Default is all features in all
            files

        subtaskName : str, optional
            The name of the subtask

        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # call the constructor from the base class (AnalysisTask)
        super(ComputeRegionMasksSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=[])

        self.geojsonFileName = geojsonFileName
        self.outFileSuffix = outFileSuffix
        self.featureList = featureList

        parentTask.add_subtask(self)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

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

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ComputeRegionMasksSubtask, self).setup_and_check()

        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS restart file found: need at least one '
                          'restart file to perform remapping of '
                          'climatologies.')

        maskSubdirectory = build_config_full_path(self.config, 'output',
                                                  'maskSubdirectory')
        make_directories(maskSubdirectory)

        if self.featureList is None:
            # get a list of features for use by other tasks (e.g. to determine
            # plot names)
            self.featureList = get_feature_list(self.geojsonFileName)

        mpasMeshName = self.config.get('input', 'mpasMeshName')

        # first, see if we have cached a mask file name in the region masks
        # directory
        regionMaskDirectory = build_config_full_path(self.config,
                                                     'diagnostics',
                                                     'regionMaskSubdirectory')
        self.maskFileName = '{}/{}_{}.nc'.format(regionMaskDirectory,
                                                 mpasMeshName,
                                                 self.outFileSuffix)

        if not os.path.exists(self.maskFileName):
            # no cached mask file, so let's see if there's already one in the
            # masks subfolder of the output directory

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            self.maskFileName = '{}/{}_{}.nc'.format(maskSubdirectory,
                                                     mpasMeshName,
                                                     self.outFileSuffix)

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        compute_region_masks(self.geojsonFileName, self.restartFileName,
                             self.maskFileName, self.featureList, self.logger)

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
