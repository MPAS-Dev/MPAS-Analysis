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

from mpas_analysis.shared.analysis_task import AnalysisTask
from mpas_analysis.shared.regions.compute_region_masks_subtask \
    import ComputeRegionMasksSubtask


class ComputeRegionMasks(AnalysisTask):
    """
    An analysis tasks for computing cell masks for regions defined by geojson
    features

    Attributes
    ----------
    regionMaskSubtasks : dict of ``ComputeRegionMasksSubtask`` objects
        The subtasks of this task with file names as keys
    """

    def __init__(self, config, conponentName):
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        conponentName : str
            The component to make mapping files for
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionMasks, self).__init__(
            config=config,
            taskName='computeRegionMasks',
            componentName=conponentName,
            tags=[])

        self.regionMaskSubtasks = {}

    def add_mask_subtask(self, regionGroup, obsFileName=None, lonVar='lon',
                         latVar='lat', meshName=None, useMpasMaskCreator=False):
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        regionGroup : str
            The name of one of the supported region groups (see
            :py:func:`mpas_analysis.shared.regions.get_region_mask()`)

        obsFileName : str, optional
            The name of an observations file to create masks for.  But default,
            lon/lat are taken from an MPAS restart file

        lonVar, latVar : str, optional
            The name of the longitude and latitude variables in ``obsFileName``

        meshName : str, optional
            The name of the mesh or grid, used as part of the mask file name.
            Default is the MPAS mesh name

        useMpasMaskCreator : bool, optional
            If ``True``, the mask creator from ``mpas_tools`` will be used
            to create the mask.  Otherwise, python code is used.  Since
            masks for observations can only be produced with the python code,
            this option is ignored if obsFileName is not ``None``.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        config = self.config

        if meshName is None:
            meshName = config.get('input', 'mpasMeshName')

        key = '{} {}'.format(meshName, regionGroup)

        if key not in self.regionMaskSubtasks:

            subprocessCount = config.getWithDefault('execute',
                                                    'parallelTaskCount',
                                                    default=1)

            if obsFileName is not None:
                useMpasMaskCreator = False

            if useMpasMaskCreator:
                subprocessCount = 1

            maskSubtask = ComputeRegionMasksSubtask(
                self, regionGroup=regionGroup, meshName=meshName,
                subprocessCount=subprocessCount, obsFileName=obsFileName,
                lonVar=lonVar, latVar=latVar,
                useMpasMaskCreator=useMpasMaskCreator)

            self.add_subtask(maskSubtask)

            self.regionMaskSubtasks[key] = maskSubtask

        return self.regionMaskSubtasks[key]
