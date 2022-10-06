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
from pyremap import LatLonGridDescriptor

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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime


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
                   'titleName': 'Significant Wave Height',
                   'era5' : 'swh'},
                  {'prefix': 'peakWaveFrequency',
                   'mpas': 'timeMonthly_avg_peakWaveFrequency',
                   'units': r's',
                   'titleName': 'Peak Wave Period',
                   'era5' : 'pp1d'},
                  {'prefix': 'iceFraction',
                   'mpas': 'timeMonthly_avg_iceFraction',
                   'units': r'',
                   'titleName': 'Ice Fraction'},
                   ]

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
        fieldsRequested = [field for field in fields if field['prefix'] in fieldList]

        variableListMpas = [field['mpas'] for field in fieldsRequested] 
        variableListObs  = [field['era5'] for field in fieldsRequested] 

        climStartYear = config.getint(sectionName, 'obsStartYear')
        climEndYear = config.getint(sectionName, 'obsEndYear')

        makeTables = config.getboolean(sectionName, 'makeTables')
        

        if makeTables:
            #for regionGroup in regionGroups:
            for season in seasons:
                tableSubtask = WavesTableSubtask(
                    parentTask=self,
                    mpasClimatologyTask=mpasClimatologyTask,
                    controlConfig=controlConfig,
                    regionMasksTask=regionMasksTask,
                    fields=fields,
                    season=season) 
                self.add_subtask(tableSubtask)


        # the variableList variables  will be added to
        # mpasClimatologyTask along with the seasons.
        #remapClimatologySubtask = RemapMpasClimatologySubtask(
        remapClimatologySubtask = RemapMpasWavesClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='wave',
            variableList=variableListMpas,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=None)

        if controlConfig is None:
            observationsDirectory = build_obs_path(
                config, 'ocean', '{}Subdirectory'.format('wave'))
            obsFileName = "{}/ERA5_1950_1978_new.nc".format(observationsDirectory)

            remapObservationsSubtask = RemapObservedWaveClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='wave',
                comparisonGridNames=comparisonGridNames,
                variableList=variableListObs)

        else:
            remapObservationsSubtask = None



        for field in fieldsRequested:

            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            configSectionName = '{}{}'.format(self.taskName, upperFieldPrefix)

            if controlConfig is None:

                refTitleLabel = "Observations (ERA 5) {:04d}-{:04d}".format(
                    climStartYear, climEndYear)
                refFieldName = field['era5']
                diffTitleLabel = None 

            else:
                controlRunName = controlConfig.get('runs', 'mainRunName')
                refTitleLabel = f'{field["titleName"]} (Control: {controlRunName})'
                diffTitleLabel = 'Main - Control'
                refFieldName = field['mpas']

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
                        remapObsClimatologySubtask=remapObservationsSubtask,
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


class RemapMpasWavesClimatology(RemapMpasClimatologySubtask): # {{{


    def customize_masked_climatology(self, climatology, season):  # {{{

        climatology['timeMonthly_avg_peakWaveFrequency'] = 1.0/climatology['timeMonthly_avg_peakWaveFrequency']

        return climatology

        # }}}
    # }}}


class RemapObservedWaveClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping waves observations
    """
    # Authors
    # -------
    # Steven Brus, Xylar Asay-Davis
    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 variableList,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        variableList : list
            List of observational variables to remap

        comparisonGridNames : list of {'latlon', 'antarctic', 'arctic'},
            optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask
        '''
        # Authors
        # -------
        # Steven Brus

        self.variableList = variableList

        # call the constructor from the base class (AnalysisTask)
        super().__init__(parentTask, seasons, fileName, outFilePrefix,
                         comparisonGridNames=['latlon'],
                         subtaskName='remapObservations')
        # }}}

    def get_observation_descriptor(self, fileName):  # {{{
        '''
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        '''
        # Authors
        # -------
        # Steven Brus 

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')
        return obsDescriptor  # }}}

    def build_observational_dataset(self, fileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Parameters
        ----------
        fileName : str
            observation file name

        Returns
        -------
        dsObs : ``xarray.Dataset``
            The observational dataset
        '''
        # Authors
        # -------
        # Steven Brus

        sectionName = self.taskName
        climStartYear = self.config.getint(sectionName, 'obsStartYear')
        climEndYear = self.config.getint(sectionName, 'obsEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_dataset(fileName, decode_times=False)
        units, reference_date = dsObs.time.attrs['units'].split('since')
        dsObs['time'] = pd.date_range(start=reference_date, periods=dsObs.sizes['time'], freq='MS')
        dsObs = dsObs.rename({'time': 'Time'})
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        dsObs = dsObs[self.variableList]

        return dsObs  # }}}

    # }}}


class WavesTableSubtask(AnalysisTask):
    def __init__(self, parentTask, mpasClimatologyTask, controlConfig,
                 regionMasksTask, season, fields, subtaskName=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ClimatologyMapWaves
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        mpasClimatologyTask : MpasClimatologyTask
            The task that produced the climatology to be remapped and plotted

        controlConfig :  MpasAnalysisConfigParser
            Configuration options for a control run (if any)

        regionMasksTask : ComputeRegionMasks
            A task for computing region masks

        fields : list of dict
            Field information 

        season : str
            One of the seasons in ``constants.monthDictionary``

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis
        tags = ['climatology', 'table']

        if subtaskName is None:
            subtaskName = f'table{season}'

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=tags)

        config = parentTask.config
        self.season = season
        self.mpasClimatologyTask = mpasClimatologyTask
        self.controlConfig = controlConfig
        self.fields = fields

        self.masksSubtask = regionMasksTask.add_mask_subtask(
            regionGroup='Arctic Sea Ice Regions')
        #self.iceShelfMasksFile = self.masksSubtask.geojsonFileName

        self.run_after(self.masksSubtask)
        self.run_after(mpasClimatologyTask)
        # }}}

    def run_task(self):  # {{{
        """
        Computes and plots table of wave climatologies 
        """
        # Authors
        # -------
        # Steven Brus
        # Xylar Asay-Davis

        self.logger.info("Computing wave climatology table...")
        config = self.config

        sectionName = self.taskName
        #regionsInTable = config.getExpression(sectionName,
        #                                      'regionsInTable')
        #if len(regionsInTable) == 0:
        #    return

        #regionsInTable = self.masksSubtask.expand_region_names(
        #    regionsInTable)

        regionsInTable = ['Baffin Bay NSIDC',
                          'Greenland Sea',
                          'Bering Sea',
                          'Barents Sea',
                          'Beaufort Sea NSIDC',
                          'Canadian Archipelago NSIDC',
                          'Central Arctic NSIDC',
                          'Chukchi Sea NSIDC',
                          'East Siberian Sea NSIDC',
                          'Hudson Bay NSIDC',
                          'Kara Sea',
                          'Laptev Sea NSIDC']

        waveFileName = get_masked_mpas_climatology_file_name(
            config, self.season, self.componentName,
            climatologyName='waveClimatologyTable')

        if not os.path.exists(waveFileName):
            with dask.config.set(schedular='threads',
                                 pool=ThreadPool(1)):

                # Load data:
                inFileName = self.mpasClimatologyTask.get_file_name(self.season)
                dsIn = xr.open_dataset(inFileName)

                ds = xr.Dataset()
                for field in self.fields:
                    mpasFieldName = field['mpas']
                    var = dsIn[mpasFieldName]
                    if 'Time' in var.dims:
                        var = var.isel(Time=0)

                    regionMaskFileName = self.masksSubtask.maskFileName
                    dsRegionMask = xr.open_dataset(regionMaskFileName)

                    # figure out the indices of the regions to plot
                    regionNames = decode_strings(dsRegionMask.regionNames)

                    regionIndices = []
                    for region in regionsInTable:
                        for index, regionName in enumerate(regionNames):
                            if region == regionName:
                                regionIndices.append(index)
                                break

                    # select only those regions we want to plot
                    dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
                    cellMasks = dsRegionMask.regionCellMasks.chunk({'nRegions': 10})

                    restartFileName = \
                        self.runStreams.readpath('restart')[0]

                    dsRestart = xr.open_dataset(restartFileName)
                    areaCell = dsRestart.areaCell

                    regionAverage = (cellMasks * areaCell * var).sum(dim='nCells')
                    regionAverage.compute()

                    totalArea = (cellMasks * areaCell).sum(dim='nCells')
                    regionAverage = regionAverage/totalArea
                    regionAverage.compute()

                    ds[mpasFieldName] = regionAverage
                    ds[mpasFieldName].attrs['units'] = field['units']
                    ds[mpasFieldName].attrs['description'] = field['titleName']

                ds['area'] = 1e-6*totalArea
                ds.area.attrs['units'] = 'km$^2$'
                ds.area.attrs['description'] = \
                    'Region area'


                ds['regionNames'] = dsRegionMask.regionNames

                write_netcdf(ds, waveFileName)
        else:
            ds = xr.open_dataset(waveFileName)

        mainRunName = config.get('runs', 'mainRunName')
        fieldNames = ['Region', 'Area', mainRunName]

        controlConfig = self.controlConfig
        if controlConfig is not None:
            controlFileName = get_masked_mpas_climatology_file_name(
                controlConfig, self.season, self.componentName,
                climatologyName='wavesTable')
            dsControl = xr.open_dataset(controlFileName)
            controlRunName = controlConfig.get('runs', 'mainRunName')
            fieldNames.append(controlRunName)
        else:
            dsControl = None
            controlRunName = None

        regionNames = decode_strings(ds.regionNames)

        outDirectory = '{}/waves/'.format(
            build_config_full_path(config, 'output', 'tablesSubdirectory'))

        try:
            os.makedirs(outDirectory)
        except OSError:
            pass

        for field in self.fields:

            tableFileName = '{}/wavesTable_{}_{}.csv'.format(outDirectory,
                                                             field['prefix'],
                                                             self.season)

            with open(tableFileName, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=fieldNames)
  
                writer.writeheader()
                for index, regionName in enumerate(regionNames):
                    if field['mpas'] == 'timeMonthly_avg_peakWaveFrequency':
                        var = 1.0/ds[field['mpas']][index].values
                    else:
                        var = ds[field['mpas']][index].values
                    row = {'Region': regionName,
                           'Area': '{}'.format(ds.area[index].values),
                           mainRunName: '{}'.format(var)}
                    if dsControl is not None:
                        row[controlRunName] = \
                            '{}'.format(dsControl[field['mpas']][index].values)
                    writer.writerow(row)

        ##### Histograms ##### 

        baseDirectory = config.get('input','baseDirectory')
        runSubdirectory = config.get('input','runSubdirectory')
        startYear = int(config.get('climatology','startYear'))
        endYear = int(config.get('climatology','endYear'))

        seasonMonths = {'ANN':[1,2,3,4,5,6,7,8,9,10,11,12], 
                        'JFM':[1,2,3],
                        'JAS':[7,8,9]}

        months = seasonMonths[self.season]

        filePrefix = baseDirectory+'/'+runSubdirectory+'/'+mainRunName+'.mpaso.waveOutput.'
        files = []
        for year in range(startYear,endYear+1):
            for month in months:
                files.append(filePrefix+str(year)+'-'+"{:02d}".format(month)+'-01_00.00.00.nc')
          

        ds = xr.open_mfdataset(files,
                               chunks={'Time': 100},
                               combine='nested',
                               concat_dim='Time',
                               parallel=True)

        regionMaskFileName = self.masksSubtask.maskFileName
        dsRegionMask = xr.open_dataset(regionMaskFileName)
   
        # figure out the indices of the regions to plot
        regionNames = decode_strings(dsRegionMask.regionNames)
   
        regionIndices = []
        for region in regionsInTable:
            for index, regionName in enumerate(regionNames):
                if region == regionName:
                    regionIndices.append(index)
                    break
   
        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
        cellMasks = dsRegionMask.regionCellMasks.chunk({'nRegions': 12})
        cellMasks = cellMasks.where(cellMasks > 0)

        for field in self.fields:
            print('computing historgram: '+field['prefix']+'_'+self.season)
            with dask.config.set(schedular='threads',
                                 pool=ThreadPool(1)):
                if field['mpas'] == 'timeMonthly_avg_peakWaveFrequency':
                    maskedVar = (cellMasks*(1.0/ds[field['prefix']])).compute()
                    varRange=(0.0,15.0)
                else: 
                    maskedVar = (cellMasks*ds[field['prefix']]).compute()
                    varRange=(0.0,10.0)

            for index, regionName in enumerate(regionNames):
                print('   '+regionName)
    
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)        
                values = maskedVar.isel(nRegions=index).values.ravel()
                ax.hist(values,range=varRange,bins=20,density=True)
                ax.set_title(regionName+' mean:'+'{:4.2f}'.format(np.nanmean(values)))
                ax.set_ylabel('Probability density')
                ax.set_xlabel(field['units'])
                plt.savefig(outDirectory+'/'+field['prefix']+'_'+self.season+'_'+regionName.replace(' ','_')+'.png')
                f = open(outDirectory+'/'+field['prefix']+'_'+self.season+'_'+regionName.replace(' ','_')+'.txt','w')
                heights = [patch.get_height() for patch in ax.patches]
                widths = [patch.get_width() for patch in ax.patches]
                xloc = [patch.get_xy()[0] for patch in ax.patches]
                for h,w,x in zip(heights,widths,xloc):
                   f.write('{} {} {}\n'.format(x,w,h))
                f.close()
                plt.close()
        ds.close()        



 # vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
