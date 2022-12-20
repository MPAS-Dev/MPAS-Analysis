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

import xarray as xr
import numpy
import matplotlib.pyplot as plt

from mpas_tools.cime.constants import constants as cime_constants

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.ocean.compute_anomaly_subtask import ComputeAnomalySubtask
from mpas_analysis.ocean.plot_hovmoller_subtask import PlotHovmollerSubtask
from mpas_analysis.ocean.plot_depth_integrated_time_series_subtask import \
    PlotDepthIntegratedTimeSeriesSubtask

from mpas_analysis.shared.constants import constants as mpas_constants


class TimeSeriesOHCAnomaly(AnalysisTask):
    """
    Performs analysis of ocean heat content (OHC) from time-series output.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    def __init__(self, config, mpasTimeSeriesTask, controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        controlconfig : mpas_tools.config.MpasConfigParser, optional
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
            tags=['timeSeries', 'ohc', 'publicObs', 'anomaly'])

        sectionName = 'timeSeriesOHCAnomaly'
        regionNames = config.getexpression(sectionName, 'regions')
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
            plotTask = PlotOHCAnomaly(
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

    def _compute_ohc(self, ds):
        """
        Compute the OHC time series.
        """

        # regionNames = self.config.getexpression('regions', 'regions')
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

        return ds


class PlotOHCAnomaly(PlotDepthIntegratedTimeSeriesSubtask):
    def customize_fig(self, fig):
        """
        A function to override to customize the figure.

        fig : matplotlib.pyplot.Figure
            The figure
        """
        def joules_to_watts_m2(joules):
            watts_m2 = joules/factor
            return watts_m2

        def watts_m2_to_joules(watts_m2):
            joules = factor*watts_m2
            return joules

        # add an axis on the right-hand side
        color = 'tab:blue'
        ax = plt.gca()
        xlim = ax.get_xlim()

        earth_surface_area = (4. * numpy.pi *
                              cime_constants['SHR_CONST_REARTH']**2)

        max_time = xlim[-1]*mpas_constants.sec_per_day

        factor = earth_surface_area*max_time/10**22

        secaxy = ax.secondary_yaxis(
            'right', functions=(joules_to_watts_m2, watts_m2_to_joules))
        secaxy.set_ylabel(r'W/m$^2$', color=color)
        secaxy.tick_params(axis='y', colors=color)
        ax.spines['right'].set_color(color)
        plt.draw()
        yticks = secaxy.get_yticks()
        for ytick in yticks:
            plt.plot(xlim, [0, watts_m2_to_joules(ytick)], color=color,
                     linewidth=0.5)
