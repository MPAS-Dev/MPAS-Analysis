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

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.mpas_xarray import mpas_xarray


def get_feature_list(config, geojsonFileName):
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
            self.featureList = get_feature_list(self.config,
                                                self.geojsonFileName)

        mpasMeshName = self.config.get('input', 'mpasMeshName')

        # first, see if we have cached a mask file name in the region masks
        # directory
        regionMaskDirectory = build_config_full_path(self.config,
                                                     'diagnostics',
                                                     'regionMaskSubdirectory')
        self.maskFileName = '{}/{}_{}.nc'.format(regionMaskDirectory,
                                                 mpasMeshName,
                                                 self.outFileSuffix)

        if os.path.exists(self.maskFileName):
            self.maskExists = True
        else:
            # no cached mask file, so let's see if there's already one in the
            # masks subfolder of the output directory

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            self.maskFileName = '{}/{}_{}.nc'.format(maskSubdirectory,
                                                     mpasMeshName,
                                                     self.outFileSuffix)

            self.maskExists = os.path.exists(self.maskFileName)

        # }}}

    def run_task(self):  # {{{
        '''
        Compute the requested climatologies
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.maskExists:
            return

        self.logger.info('Creating masks file {}'.format(self.maskFileName))

        with xr.open_dataset(self.restartFileName) as dsRestart:
            dsRestart = mpas_xarray.subset_variables(dsRestart,
                                                     ['lonCell', 'latCell'])
            latCell = numpy.rad2deg(dsRestart.latCell.values)

            # transform longitudes to [-180, 180)
            lonCell = numpy.mod(numpy.rad2deg(dsRestart.lonCell.values) + 180.,
                                360.) - 180.

        # create shapely geometry for lonCell and latCell
        cellPoints = [shapely.geometry.Point(x, y) for x, y in
                      zip(lonCell, latCell)]

        nCells = len(cellPoints)

        masks = []
        regionNames = []
        nChar = 0
        self.logger.info('  Computing masks from {}...'.format(
                self.geojsonFileName))
        with open(self.geojsonFileName) as f:
            featureData = json.load(f)

        for feature in featureData['features']:
            name = feature['properties']['name']
            if name not in self.featureList:
                continue

            self.logger.info('      {}'.format(name))

            shape = shapely.geometry.shape(feature['geometry'])
            mask = numpy.array([shape.contains(point) for point
                                in cellPoints], dtype=bool)

            nChar = max(nChar, len(name))

            masks.append(mask)
            regionNames.append(name)

        # create a new data array for masks and another for mask names
        self.logger.info('  Creating and writing masks dataset...')
        nRegions = len(regionNames)
        dsMasks = xr.Dataset()
        dsMasks['masks'] = (('nRegions', 'nCells'),
                            numpy.zeros((nRegions, nCells), dtype=bool))
        dsMasks['regionNames'] = (('nRegions'),
                                  numpy.zeros((nRegions),
                                              dtype='|S{}'.format(nChar)))
        for index in range(nRegions):
            regionName = regionNames[index]
            mask = masks[index]
            dsMasks['regionCellMasks'][index, :] = mask
            dsMasks['regionNames'][index] = regionName

        write_netcdf(dsMasks, self.maskFileName)

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
