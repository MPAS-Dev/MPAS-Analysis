# -*- coding: utf-8 -*-
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
#

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.compute_anomaly_subtask import ComputeAnomalySubtask
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask

from mpas_analysis.ocean.utility import get_standard_region_names


class TimeSeriesSalinityAnomaly(AnalysisTask):
    """
    Performs analysis of time series of salinity anomalies from the first
    simulation year as a function of depth.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasTimeSeriesTask):
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  mpas_tools.config.MpasConfigParser
            Contains configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesSalinityAnomaly, self).__init__(
            config=config,
            taskName='timeSeriesSalinityAnomaly',
            componentName='ocean',
            tags=['timeSeries', 'salinity', 'publicObs', 'anomaly'])

        sectionName = 'hovmollerSalinityAnomaly'
        regionShortNames = config.getexpression(sectionName,
                                                'regionShortNames')
        movingAveragePoints = config.getint(sectionName, 'movingAveragePoints')

        mpasFieldName = 'timeMonthly_avg_avgValueWithinOceanLayerRegion_' \
            'avgLayerSalinity'

        timeSeriesFileName = 'regionAveragedSalinityAnomaly.nc'

        anomalyTask = ComputeAnomalySubtask(
            parentTask=self,
            mpasTimeSeriesTask=mpasTimeSeriesTask,
            outFileName=timeSeriesFileName,
            variableList=[mpasFieldName],
            movingAveragePoints=movingAveragePoints)
        self.add_subtask(anomalyTask)

        regionNames = get_standard_region_names(config, regionShortNames)

        for regionName in regionNames:
            caption = f'Trend of {regionName} Salinity Anomaly vs depth'
            plotTask = PlotHovmollerSubtask(
                parentTask=self,
                regionName=regionName,
                inFileName=timeSeriesFileName,
                outFileLabel='SAnomalyZ',
                fieldNameInTitle='Salinity Anomaly',
                mpasFieldName=mpasFieldName,
                unitsLabel='[PSU]',
                sectionName=sectionName,
                thumbnailSuffix=u'ΔS',
                imageCaption=caption,
                galleryGroup='Trends vs Depth',
                groupSubtitle=None,
                groupLink='trendsvsdepth',
                galleryName=None)

            plotTask.run_after(anomalyTask)
            self.add_subtask(plotTask)
