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
import csv
import xarray as xr
import dask
from multiprocessing.pool import ThreadPool
from pyremap import ProjectionGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path, decode_strings, \
    build_config_full_path
from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask, get_antarctic_stereographic_projection
from mpas_analysis.shared.climatology.climatology import \
    get_masked_mpas_climatology_file_name

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask    

from mpas_analysis.shared.constants import constants


class ClimatologyMapWaves(AnalysisTask):  # {{{
    """
    An analysis task for the computation of wave climatologies
    """
    # Authors
    # -------
    # Steven Brus
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask,
                 controlConfig):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig :  ``MpasAnalysisConfigParser``
            Configuration options for a control run
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis

        fields = [{'prefix': 'significantWaveHeight',
                   'mpas': 'timeMonthly_avg_significantWaveHeight',
                   'units': r'm',
                   'titleName': 'Significant Wave Height'},
                  {'prefix': 'peakWaveFrequency',
                   'mpas': 'timeMonthly_avg_peakWaveFrequency',
                   'units': r's$^{-1}$',
                   'titleName': 'Peak Wave Frequency'}]

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config, taskName='climatologyMapWaves',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'waves'])

        sectionName = self.taskName

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        regionGroups = config.getExpression(sectionName, 'regionGroups')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        fieldList = config.getExpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        variableList = [field['mpas'] for field in fields] 

        makeTables = config.getboolean(sectionName, 'makeTables')
        
        print(variableList)

        #if makeTables:
        #    for regionGroup in regionGroups:
        #        for season in seasons:
        #            tableSubtask = WavesTableSubtask(
        #                parentTask=self,
        #                mpasClimatologyTask=mpasClimatologyTask,
        #                controlConfig=controlConfig,
        #                regionMasksTask=regionMasksTask,
        #                regionGroup=regionGroup,
        #                fields=fields,
        #                season=season) 
        #            self.add_subtask(tableSubtask)


        # the variableList variables  will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='wave',
            variableList=variableList,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=None)


        for field in fields:

            print(field)
            fieldPrefix = field['prefix']
            refFieldName = field['mpas']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            configSectionName = '{}{}'.format(self.taskName, upperFieldPrefix)

            if controlConfig is None:

                galleryName = field['titleName']
                refTitleLabel = None 
                diffTitleLabel = None 

            else:
                controlRunName = controlConfig.get('runs', 'mainRunName')
                refTitleLabel = f'{field["titleName"]} (Control: {controlRunName})'
                diffTitleLabel = 'Main - Control'

            outFileLabel = '{}Wave'.format(fieldPrefix)

            for comparisonGridName in comparisonGridNames:
                for season in seasons:

                    subtaskName = 'plot{}_{}_{}'.format(upperFieldPrefix,
                                                        season,
                                                        comparisonGridName)

                    # make a new subtask for this season and comparison grid
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self, 
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        controlConfig=controlConfig,
                        subtaskName=subtaskName)


                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=field['titleName'],
                        mpasFieldName=field['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=field['units'],
                        imageCaption=field['titleName'],
                        #galleryGroup=field['titleName'],
                        galleryGroup='waves',
                        groupSubtitle=None,
                        #groupLink='{}Waves'.format(fieldPrefix),
                        groupLink='waves',
                        galleryName=field['titleName'],
                        configSectionName=configSectionName)

                    self.add_subtask(subtask)
        # }}}

    def setup_and_check(self):  # {{{
       """
       Perform steps to set up the analysis and check for errors in the setup.
       """
       # Authors
       # -------
       # Steven Brus
       # Xylar Asay-Davis

       # first, call setup_and_check from the base class
       # (AnalysisTask), which will perform some common setup
       super().setup_and_check()

       if not self.namelist.getbool('config_use_active_wave'):
        
           raise ValueError('*** climatologyMapWaves requires '
                            'config_use_active_wave\n'
                            '    to be .true.  Otherwise, no '
                            'wave data is available \n'
                            '    for plotting.')
        # }}}
    # }}}




#class WavesTableSubtask(AnalysisTask):
#    def __init__(self, parentTask, mpasClimatologyTask, controlConfig,
#                 regionMasksTask, regionGroup, season, fields, subtaskName=None):  # {{{
#        """
#        Construct the analysis task.
#
#        Parameters
#        ----------
#        parentTask :  ClimatologyMapWaves
#            The parent task, used to get the ``taskName``, ``config`` and
#            ``componentName``
#
#        mpasClimatologyTask : MpasClimatologyTask
#            The task that produced the climatology to be remapped and plotted
#
#        controlConfig :  MpasAnalysisConfigParser
#            Configuration options for a control run (if any)
#
#        regionMasksTask : ComputeRegionMasks
#            A task for computing region masks
#
#        regionGroup : str 
#            Region group name where each region is a row of the table
#
#        fields : list of dict
#            Field information 
#
#        season : str
#            One of the seasons in ``constants.monthDictionary``
#
#        subtaskName : str, optional
#            The name of the subtask
#        """
#        # Authors
#        # -------
#        # Steven Brus
#        # Xylar Asay-Davis
#        tags = ['climatology', 'table']
#
#        if subtaskName is None:
#            subtaskName = f'table{season}{regionGroup.replace(" ","")}'
#
#        # call the constructor from the base class (AnalysisTask)
#        super().__init__(
#            config=parentTask.config,
#            taskName=parentTask.taskName,
#            subtaskName=subtaskName,
#            componentName=parentTask.componentName,
#            tags=tags)
#
#        config = parentTask.config
#        self.season = season
#        self.mpasClimatologyTask = mpasClimatologyTask
#        self.controlConfig = controlConfig
#        self.fields = fields
#
#        self.masksSubtask = regionMasksTask.add_mask_subtask(
#             regionGroup=regionGroup)
#        self.iceShelfMasksFile = self.masksSubtask.geojsonFileName
#
#        self.run_after(self.masksSubtask)
#        self.run_after(mpasClimatologyTask)
#        # }}}
#
#    def run_task(self):  # {{{
#        """
#        Computes and plots table of wave climatologies 
#        """
#        # Authors
#        # -------
#        # Steven Brus
#        # Xylar Asay-Davis
#
#        self.logger.info("Computing wave climatology table...")
#        config = self.config
#
#        sectionName = self.taskName
#        regionsInTable = config.getExpression(sectionName,
#                                                 'regionsInTable')
#        if len(regionsInTable) == 0:
#            return
#
#        regionsInTable = self.masksSubtask.expand_region_names(
#            regionsInTable)
#
#        waveFileName = get_masked_mpas_climatology_file_name(
#            config, self.season, self.componentName,
#            climatologyName='waveClimatologyTable')
#
#        if not os.path.exists(waveFileName):
#            with dask.config.set(schedular='threads',
#                                 pool=ThreadPool(1)):
#
#                # Load data:
#                inFileName = self.mpasClimatologyTask.get_file_name(self.season)
#                dsIn = xr.open_dataset(inFileName)
#
#                ds = xr.Dataset()
#                for field in self.fields:
#                    mpasFieldName = field['mpas']
#                    var = dsIn[mpasFieldName]
#                    if 'Time' in var.dims:
#                        var = var.isel(Time=0)
#
#                    regionMaskFileName = self.masksSubtask.maskFileName
#                    dsRegionMask = xr.open_dataset(regionMaskFileName)
#
#                    # figure out the indices of the regions to plot
#                    regionNames = decode_strings(dsRegionMask.regionNames)
#
#                    regionIndices = []
#                    for region in regionsInTable:
#                        for index, regionName in enumerate(regionNames):
#                            if region == regionName:
#                                regionIndices.append(index)
#                                break
#
#                    # select only those regions we want to plot
#                    dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
#                    cellMasks = dsRegionMask.regionCellMasks.chunk({'nRegions': 10})
#
#                    restartFileName = \
#                        self.runStreams.readpath('restart')[0]
#
#                    dsRestart = xr.open_dataset(restartFileName)
#                    areaCell = dsRestart.areaCell
#
#                    regionAverage = (areaCell * field).sum(dim='nCells')
#                    regionAverage.compute()
#                    regionAverage = regionAverage/totalArea
#
#                    ds[mpasFieldName] = var
#                    ds[mpasFieldName].attrs['units'] = field['units']
#                    ds[mpasFieldName].attrs['description'] = field['titleName']
#
#               ds['area'] = 1e-6*totalArea
#               ds.area.attrs['units'] = 'km$^2$'
#               ds.area.attrs['description'] = \
#                   'Region or ice shelf area'
#
#
#               ds['regionNames'] = dsRegionMask.regionNames
#
#               write_netcdf(ds, meltRateFileName)
#        else:
#            ds = xr.open_dataset(meltRateFileName)
#
#        mainRunName = config.get('runs', 'mainRunName')
#        fieldNames = ['Region', 'Area', mainRunName]
#
#        controlConfig = self.controlConfig
#        if controlConfig is not None:
#            controlFileName = get_masked_mpas_climatology_file_name(
#                controlConfig, self.season, self.componentName,
#                climatologyName='antarcticMeltTable')
#            dsControl = xr.open_dataset(controlFileName)
#            controlRunName = controlConfig.get('runs', 'mainRunName')
#            fieldNames.append(controlRunName)
#        else:
#            dsControl = None
#            controlRunName = None
#
#        regionNames = decode_strings(ds.regionNames)
#
#        outDirectory = '{}/antarcticMelt/'.format(
#            build_config_full_path(config, 'output', 'tablesSubdirectory'))
#
#        try:
#            os.makedirs(outDirectory)
#        except OSError:
#            pass
#
#        tableFileName = '{}/antarcticMeltRateTable_{}.csv'.format(outDirectory,
#                                                                  self.season)
#
#        with open(tableFileName, 'w', newline='') as csvfile:
#            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)
#
#            writer.writeheader()
#            for index, regionName in enumerate(regionNames):
#                row = {'Region': regionName,
#                       'Area': '{}'.format(ds.area[index].values),
#                       mainRunName: '{}'.format(ds.meltRates[index].values)}
#                if dsControl is not None:
#                    row[controlRunName] = \
#                        '{}'.format(dsControl.meltRates[index].values)
#                writer.writerow(row)
#
#        tableFileName = '{}/antarcticMeltFluxTable_{}.csv'.format(outDirectory,
#                                                                  self.season)
#
#        with open(tableFileName, 'w', newline='') as csvfile:
#            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)
#
#            writer.writeheader()
#            for index, regionName in enumerate(regionNames):
#                row = {'Region': regionName,
#                       'Area': '{}'.format(ds.area[index].values),
#                       mainRunName: '{}'.format(ds.totalMeltFlux[index].values)}
#                if dsControl is not None:
#                    row[controlRunName] = \
#                        '{}'.format(dsControl.totalMeltFlux[index].values)
#                writer.writerow(row)
#
#        # }}}
#    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
