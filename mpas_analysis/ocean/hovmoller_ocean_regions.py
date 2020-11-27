# -*- coding: utf-8 -*-
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
#

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask


class HovmollerOceanRegions(AnalysisTask):  # {{{
    """
    Compute and plot a Hovmoller diagram (depth vs. time) for regionally
    analyzed data.  The mean of the data are computed over each region at each
    depth and time.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, regionMasksTask, oceanRegionalProfilesTask,
                 controlConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        oceanRegionalProfilesTask : mpas_analysis.ocean.OceanRegionalProfiles
            A task for computing ocean regional profiles

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(HovmollerOceanRegions, self).__init__(
            config=config,
            taskName='hovmollerOceanRegions',
            componentName='ocean',
            tags=['profiles', 'timeSeries', 'hovmoller'])

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.get('timeSeries', 'endYear')
        if endYear == 'end':
            # a valid end year wasn't found, so likely the run was not found,
            # perhaps because we're just listing analysis tasks
            endYear = startYear
        else:
            endYear = int(endYear)

        regionGroups = config.getExpression('hovmollerOceanRegions',
                                            'regionGroups')

        for regionGroup in regionGroups:
            suffix = regionGroup[0].upper() + regionGroup[1:].replace(' ', '')
            regionGroupSection = 'hovmoller{}'.format(suffix)
            regionNames = config.getExpression(regionGroupSection,
                                               'regionNames')
            if len(regionNames) == 0:
                return

            fields = config.getExpression(regionGroupSection, 'fields')

            masksSubtask = regionMasksTask.add_mask_subtask(regionGroup)
            masksFile = masksSubtask.geojsonFileName
            timeSeriesName = masksSubtask.outFileSuffix

            regionNames = masksSubtask.expand_region_names(regionNames)

            self.masksSubtask = masksSubtask

            oceanRegionalProfilesTask.add_region_group(
                regionMasksTask, regionGroup, regionNames,
                fields, startYear, endYear)

            combineSubtask = oceanRegionalProfilesTask.combineSubtasks[
                regionGroup][(startYear, endYear)]

            for field in fields:
                prefix = field['prefix']
                suffix = prefix[0].upper() + prefix[1:]
                for regionName in regionNames:
                    subtaskName = 'plotHovmoller_{}_{}'.format(
                        prefix, regionName.replace(' ', '_'))
                    inFileName = \
                        '{}/regionalProfiles_{}_{:04d}-{:04d}.nc'.format(
                            timeSeriesName, timeSeriesName,
                            startYear, endYear)
                    titleName = field['titleName']
                    caption = 'Time series of {} {} vs ' \
                              'depth'.format(regionName.replace('_', ' '),
                                             titleName)
                    hovmollerSubtask = PlotHovmollerSubtask(
                        parentTask=self,
                        regionName=regionName,
                        inFileName=inFileName,
                        outFileLabel='{}_hovmoller'.format(prefix),
                        fieldNameInTitle=titleName,
                        mpasFieldName='{}_mean'.format(prefix),
                        unitsLabel=field['units'],
                        sectionName='hovmollerOceanRegions{}'.format(suffix),
                        thumbnailSuffix='',
                        imageCaption=caption,
                        galleryGroup='{} Time Series vs Depths'.format(
                            regionGroup),
                        groupSubtitle=None,
                        groupLink='ocnreghovs',
                        galleryName=titleName,
                        subtaskName=subtaskName,
                        controlConfig=controlConfig,
                        regionMaskFile=masksFile)
                    hovmollerSubtask.run_after(combineSubtask)
                    self.add_subtask(hovmollerSubtask)

        # }}}
    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
