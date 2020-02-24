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

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray as xr
from geometric_features import read_feature_collection
import mpas_tools.conversion

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, get_region_mask
from mpas_analysis.shared.io import write_netcdf


def compute_mpas_transect_masks(geojsonFileName, meshFileName, maskFileName,
                                logger=None):
    """
    Build a transect mask file from the given MPAS mesh and geojson file \
    defining a set of transects.
    """
    if os.path.exists(maskFileName):
        return

    dsMesh = xr.open_dataset(meshFileName)
    fcMask = read_feature_collection(geojsonFileName)
    dsMask = mpas_tools.conversion.mask(dsMesh=dsMesh, fcMask=fcMask,
                                        logger=logger)

    write_netcdf(dsMask, maskFileName)


class ComputeTransectMasksSubtask(AnalysisTask):  # {{{
    """
    An analysis tasks for computing cell masks for transects defined by geojson
    features

    Attributes
    ----------
    geojsonFileName : str
        A geojson file, typically from the MPAS ``geometric_features``
        repository, defining the shapes to be masked

    outFileSuffix : str
        The suffix for the resulting mask file

    maskFileName : str
        The name of the output mask file

    meshFileName : str
        A mesh file used to create the masks
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, geojsonFileName, outFileSuffix,
                 subtaskName='computeTransectMasks', subprocessCount=1):
        # {{{
        """
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

        subtaskName : str, optional
            The name of the subtask

        subprocessCount : int, optional
            The number of processes that can be used to make the mask
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # call the constructor from the base class (AnalysisTask)
        super(ComputeTransectMasksSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=[])

        self.geojsonFileName = geojsonFileName
        self.outFileSuffix = outFileSuffix
        self.subprocessCount = subprocessCount

        # }}}

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
        super(ComputeTransectMasksSubtask, self).setup_and_check()

        try:
            self.obsFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS restart file found: need at least one '
                          'restart file to perform region masking.')

        maskSubdirectory = build_config_full_path(self.config, 'output',
                                                  'maskSubdirectory')
        make_directories(maskSubdirectory)

        # first, see if we have cached a mask file name in the region masks
        # directory

        meshName = self.config.get('input', 'mpasMeshName')

        self.maskFileName = get_region_mask(
            self.config, '{}_{}.nc'.format(meshName, self.outFileSuffix))

        if not os.path.exists(self.maskFileName):
            # no cached mask file, so let's see if there's already one in the
            # masks subfolder of the output directory

            maskSubdirectory = build_config_full_path(self.config, 'output',
                                                      'maskSubdirectory')
            self.maskFileName = '{}/{}_{}.nc'.format(maskSubdirectory,
                                                     meshName,
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

        compute_mpas_transect_masks(
            self.geojsonFileName, self.obsFileName, self.maskFileName,
            logger=self.logger)

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
