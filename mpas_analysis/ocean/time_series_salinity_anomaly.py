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

from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.compute_anomaly_subtask import ComputeAnomalySubtask
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask


class TimeSeriesSalinityAnomaly(AnalysisTask):
    """
    Performs analysis of time series of salinity anomalies from the first
    simulation year as a function of depth.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasTimeSeriesTask):  # {{{
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
        regionNames = config.getExpression(sectionName, 'regions')
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

        for regionName in regionNames:
            caption = 'Trend of {} Salinity Anomaly vs depth'.format(
                regionName)
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

        # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
