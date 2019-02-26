# -*- coding: utf-8 -*-
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.compute_anomaly_subtask import ComputeAnomalySubtask
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask
from mpas_analysis.ocean.plot_depth_integrated_time_series_subtask import \
    PlotDepthIntegratedTimeSeriesSubtask


class TimeSeriesOHCAnomaly(AnalysisTask):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    def __init__(self, config, mpasTimeSeriesTask, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesOHCAnomaly, self).__init__(
            config=config,
            taskName='timeSeriesOHCAnomaly',
            componentName='ocean',
            tags=['timeSeries', 'ohc', 'publicObs'])

        sectionName = 'timeSeriesOHCAnomaly'
        regionNames = config.getExpression(sectionName, 'regions')
        movingAveragePoints = config.getint(sectionName, 'movingAveragePoints')

        self.variableDict = {}
        for suffix in ['avgLayerTemperature', 'sumLayerMaskValue',
                       'avgLayerArea', 'avgLayerThickness']:
            key = 'timeMonthly_avg_avgValueWithinOceanLayerRegion_' + suffix
            self.variableDict[key] = suffix

        mpasFieldName = 'ohc'

        timeSeriesFileName = 'regionAveragedOHCAnomaly.nc'

        anomalyTask = ComputeAnomalySubtask(
            parentTask=self,
            mpasTimeSeriesTask=mpasTimeSeriesTask,
            outFileName=timeSeriesFileName,
            variableList=list(self.variableDict.keys()),
            movingAveragePoints=movingAveragePoints,
            alter_dataset=self._compute_ohc)
        self.add_subtask(anomalyTask)

        for regionName in regionNames:
            caption = 'Trend of {} OHC Anomaly vs depth'.format(
                regionName)
            plotTask = PlotHovmollerSubtask(
                parentTask=self,
                regionName=regionName,
                inFileName=timeSeriesFileName,
                outFileLabel='ohcAnomalyZ',
                fieldNameInTitle='OHC Anomaly',
                mpasFieldName=mpasFieldName,
                unitsLabel=r'[$\times 10^{22}$ J]',
                sectionName='hovmollerOHCAnomaly',
                thumbnailSuffix=u'ΔOHC',
                imageCaption=caption,
                galleryGroup='Trends vs Depth',
                groupSubtitle=None,
                groupLink='trendsvsdepth',
                galleryName=None)

            plotTask.run_after(anomalyTask)
            self.add_subtask(plotTask)

            caption = 'Running Mean of the Anomaly in {} Ocean Heat ' \
                'Content'.format(regionName)
            plotTask = PlotDepthIntegratedTimeSeriesSubtask(
                parentTask=self,
                regionName=regionName,
                inFileName=timeSeriesFileName,
                outFileLabel='ohcAnomaly',
                fieldNameInTitle='OHC Anomaly',
                mpasFieldName=mpasFieldName,
                yAxisLabel=r'$\Delta$OHC [$\times 10^{22}$ J]',
                sectionName='timeSeriesOHCAnomaly',
                thumbnailSuffix=u'ΔOHC',
                imageCaption=caption,
                galleryGroup='Time Series',
                groupSubtitle=None,
                groupLink='timeseries',
                galleryName=None,
                controlConfig=controlConfig)

            plotTask.run_after(anomalyTask)
            self.add_subtask(plotTask)

        # }}}

    def _compute_ohc(self, ds):  # {{{
        '''
        Compute the OHC time series.
        '''

        # regionNames = self.config.getExpression('regions', 'regions')
        # ds['regionNames'] = ('nOceanRegionsTmp', regionNames)

        # for convenience, rename the variables to simpler, shorter names
        ds = ds.rename(self.variableDict)

        # specific heat [J/(kg*degC)]
        cp = self.namelist.getfloat('config_specific_heat_sea_water')
        # [kg/m3]
        rho = self.namelist.getfloat('config_density0')

        unitsScalefactor = 1e-22

        ds['ohc'] = unitsScalefactor * rho * cp * ds['sumLayerMaskValue'] * \
            ds['avgLayerArea'] * ds['avgLayerThickness'] * \
            ds['avgLayerTemperature']
        ds.ohc.attrs['units'] = '$10^{22}$ J'
        ds.ohc.attrs['description'] = 'Ocean heat content in each region'

        # Note: restart file, not a mesh file because we need refBottomDepth,
        # not in a mesh file
        try:
            restartFile = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for OHC calculation')

        # Define/read in general variables
        with xr.open_dataset(restartFile) as dsRestart:
            # reference depth [m]
            # add depths as a coordinate to the data set
            ds.coords['depth'] = (('nVertLevels',),
                                  dsRestart.refBottomDepth.values)

        return ds  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
