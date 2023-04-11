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
import os
import xarray
import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import gsw
from gsw import freezing
from gsw.density import sigma0

import dask
import multiprocessing
from multiprocessing.pool import ThreadPool

from geometric_features import FeatureCollection, read_feature_collection

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.plot import savefig, add_inset

from mpas_analysis.shared.io import write_netcdf_with_fill

from mpas_analysis.shared.io.utility import decode_strings, \
    build_obs_path, build_config_full_path, make_directories

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.ocean.utility import compute_zmid

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name, \
    get_masked_mpas_climatology_file_name

from mpas_analysis.shared.plot.colormap import register_custom_colormaps
from mpas_analysis.shared.plot.title import limit_title


class RegionalTSDiagrams(AnalysisTask):
    """
    Create T-S Diagrams of the climatology within a given ocean region

    Attributes
    ----------
    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask,
                 controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(RegionalTSDiagrams, self).__init__(
            config=config,
            taskName='regionalTSDiagrams',
            componentName='ocean',
            tags=['climatology', 'regions', 'publicObs', 'TSDiagrams'])

        self.run_after(mpasClimatologyTask)
        self.mpasClimatologyTask = mpasClimatologyTask

        regionGroups = config.getexpression(self.taskName, 'regionGroups')

        self.seasons = config.getexpression(self.taskName, 'seasons')

        obsDicts = {
            'SOSE': {
                'suffix': 'SOSE',
                'gridName': 'SouthernOcean_0.167x0.167degree',
                'gridFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                                'SouthernOcean_0.167x0.167degree_20180710.nc',
                'TFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                             'SouthernOcean_0.167x0.167degree_20180710.nc',
                'SFileName': 'SOSE/SOSE_2005-2010_monthly_salinity_'
                             'SouthernOcean_0.167x0.167degree_20180710.nc',
                'volFileName': 'SOSE/SOSE_volume_'
                               'SouthernOcean_0.167x0.167degree_20190815.nc',
                'preprocessedFileTemplate':
                    'SOSE/SOSE_{}_T_S_z_vol_'
                    'SouthernOcean_0.167x0.167degree_20200514.nc',
                'lonVar': 'lon',
                'latVar': 'lat',
                'TVar': 'theta',
                'SVar': 'salinity',
                'volVar': 'volume',
                'zVar': 'z',
                'tVar': 'Time'},
            'WOA18': {
                'suffix': 'WOA18',
                'gridName': 'Global_0.25x0.25degree',
                'gridFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'TFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'SFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
                'volFileName': None,
                'preprocessedFileTemplate':
                    'WOA18/woa18_{}_T_S_z_vol_20200514.nc',
                'lonVar': 'lon',
                'latVar': 'lat',
                'TVar': 't_an',
                'SVar': 's_an',
                'volVar': 'volume',
                'zVar': 'depth',
                'tVar': 'month'}}

        allObsUsed = []

        for regionGroup in regionGroups:
            sectionSuffix = regionGroup[0].upper() + \
                regionGroup[1:].replace(' ', '')
            sectionName = 'TSDiagramsFor{}'.format(sectionSuffix)
            obsList = config.getexpression(sectionName, 'obs')
            allObsUsed = allObsUsed + obsList
        allObsUsed = set(allObsUsed)

        for obsName in allObsUsed:
            obsDict = obsDicts[obsName]
            obsDicts[obsName]['climatologyTask'] = {}
            for season in self.seasons:

                climatologySubtask = ComputeObsTSClimatology(
                    self, obsName, obsDict, season=season)
                obsDicts[obsName]['climatologyTask'][season] = \
                    climatologySubtask
                self.add_subtask(climatologySubtask)

        for regionGroup in regionGroups:
            sectionSuffix = regionGroup[0].upper() + \
                regionGroup[1:].replace(' ', '')
            sectionName = 'TSDiagramsFor{}'.format(sectionSuffix)

            regionNames = config.getexpression(sectionName, 'regionNames')
            if len(regionNames) == 0:
                continue

            mpasMasksSubtask = regionMasksTask.add_mask_subtask(
                regionGroup=regionGroup)

            try:
                regionNames = mpasMasksSubtask.expand_region_names(regionNames)
            except FileNotFoundError:
                # this may happen if we can't create the geojson file to expand
                # its contents, e.g. if we're just doing mpas_analysis --list
                regionNames = []

            obsList = config.getexpression(sectionName, 'obs')
            groupObsDicts = {}

            for obsName in obsList:
                localObsDict = dict(obsDicts[obsName])
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                obsMasksSubtask = regionMasksTask.add_mask_subtask(
                    regionGroup, obsFileName=obsFileName,
                    lonVar=localObsDict['lonVar'],
                    latVar=localObsDict['latVar'],
                    meshName=localObsDict['gridName'])

                obsDicts[obsName]['maskTask'] = obsMasksSubtask

                localObsDict['maskTask'] = obsMasksSubtask
                groupObsDicts[obsName] = localObsDict

            for regionName in regionNames:

                fullSuffix = sectionSuffix + '_' + \
                    regionName.replace(' ', '')

                for season in self.seasons:
                    computeRegionSubtask = ComputeRegionTSSubtask(
                        self, regionGroup, regionName, controlConfig,
                        sectionName, fullSuffix, mpasClimatologyTask,
                        mpasMasksSubtask, groupObsDicts, season)
                    computeRegionSubtask.run_after(mpasMasksSubtask)
                    for obsName in obsList:
                        task = \
                            obsDicts[obsName]['climatologyTask'][season]
                        computeRegionSubtask.run_after(task)
                        task = obsDicts[obsName]['maskTask']
                        computeRegionSubtask.run_after(task)

                    self.add_subtask(computeRegionSubtask)

                    plotRegionSubtask = PlotRegionTSDiagramSubtask(
                        self, regionGroup, regionName, controlConfig,
                        sectionName, fullSuffix, mpasClimatologyTask,
                        mpasMasksSubtask, groupObsDicts, season)

                    plotRegionSubtask.run_after(computeRegionSubtask)
                    self.add_subtask(plotRegionSubtask)

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(RegionalTSDiagrams, self).setup_and_check()

        # don't add the variables and seasons to mpasClimatologyTask until
        # we're sure this subtask is supposed to run
        variableList = ['timeMonthly_avg_activeTracers_temperature',
                        'timeMonthly_avg_activeTracers_salinity',
                        'timeMonthly_avg_layerThickness']
        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=self.seasons)


class ComputeObsTSClimatology(AnalysisTask):
    """
    Compute the seasonal climatology for the given observational data set

    Attributes
    ----------
    obsDict : dicts
        Information on the observational data sets

    season : str
        The season to compute the climatology for
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, obsName, obsDict, season):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        obsDict : dicts
            Information on the observational data sets

        season : str
            The season to compute the climatology for
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeObsTSClimatology, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='climatolgy{}_{}'.format(obsDict['suffix'], season))

        self.obsName = obsName
        self.obsDict = obsDict
        self.season = season

        self.fileName = self._get_file_name(obsDict)
        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))

    def run_task(self):
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if os.path.exists(self.fileName):
            return

        config = self.config
        obsDict = self.obsDict
        season = self.season

        obsFileName = build_obs_path(
            config, component=self.componentName,
            relativePath=obsDict['preprocessedFileTemplate'].format(season))

        if os.path.exists(obsFileName):
            # we already pre-processed this data so no need to redo it!
            os.symlink(obsFileName, self.fileName)
            return

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            self.logger.info("\n computing T S climatogy for {}...".format(
                self.obsName))

            chunk = {obsDict['tVar']: 6}

            TVarName = obsDict['TVar']
            SVarName = obsDict['SVar']
            zVarName = obsDict['zVar']
            lonVarName = obsDict['lonVar']
            latVarName = obsDict['latVar']
            volVarName = obsDict['volVar']

            obsFileName = build_obs_path(
                config, component=self.componentName,
                relativePath=obsDict['TFileName'])
            self.logger.info('  Reading from {}...'.format(obsFileName))
            ds = xarray.open_dataset(obsFileName, chunks=chunk)
            if obsDict['SFileName'] != obsDict['TFileName']:
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=obsDict['SFileName'])
                self.logger.info('  Reading from {}...'.format(obsFileName))
                dsS = xarray.open_dataset(obsFileName, chunks=chunk)
                ds[SVarName] = dsS[SVarName]

            if obsDict['volFileName'] is None:
                # compute volume from lat, lon, depth bounds
                self.logger.info('  Computing volume...'.format(obsFileName))
                latBndsName = ds[latVarName].attrs['bounds']
                lonBndsName = ds[lonVarName].attrs['bounds']
                zBndsName = ds[zVarName].attrs['bounds']
                latBnds = ds[latBndsName]
                lonBnds = ds[lonBndsName]
                zBnds = ds[zBndsName]
                dLat = numpy.deg2rad(latBnds[:, 1] - latBnds[:, 0])
                dLon = numpy.deg2rad(lonBnds[:, 1] - lonBnds[:, 0])
                lat = numpy.deg2rad(ds[latVarName])
                dz = zBnds[:, 1] - zBnds[:, 0]
                radius = 6378137.0
                area = radius**2*numpy.cos(lat)*dLat*dLon
                ds[volVarName] = dz*area

            elif obsDict['volFileName'] != obsDict['TFileName']:
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=obsDict['volFileName'])
                self.logger.info('  Reading from {}...'.format(obsFileName))
                dsVol = xarray.open_dataset(obsFileName)
                ds[volVarName] = dsVol[volVarName]

            temp_file_name = self._get_file_name(obsDict, suffix='_combined')
            temp_files = [temp_file_name]
            self.logger.info('  computing the dataset')
            ds.compute()
            self.logger.info('  writing temp file {}...'.format(temp_file_name))
            write_netcdf_with_fill(ds, temp_file_name)

            chunk = {obsDict['latVar']: 400,
                     obsDict['lonVar']: 400}

            self.logger.info('  Reading back from {}...'.format(temp_file_name))
            ds = xarray.open_dataset(temp_file_name, chunks=chunk)

            if obsDict['tVar'] in ds.dims:
                if obsDict['tVar'] != 'Time':
                    if obsDict['tVar'] == 'month':
                        ds = ds.rename({obsDict['tVar']: 'Time'})
                        ds.coords['month'] = ds['Time']
                    else:
                        ds = ds.rename({obsDict['tVar']: 'Time'})
                if 'year' not in ds:
                    ds.coords['year'] = numpy.ones(ds.sizes['Time'], int)

                monthValues = constants.monthDictionary[season]
                self.logger.info('  Computing climatology...')
                ds = compute_climatology(ds, monthValues, maskVaries=True)

            temp_file_name = self._get_file_name(obsDict, suffix='_clim')
            temp_files.append(temp_file_name)
            self.logger.info('  computing the dataset')
            ds.compute()
            self.logger.info('  writing temp file {}...'.format(temp_file_name))
            write_netcdf_with_fill(ds, temp_file_name)
            self.logger.info('  Reading back from {}...'.format(temp_file_name))
            ds = xarray.open_dataset(temp_file_name, chunks=chunk)

            self.logger.info('  Broadcasting z coordinate...')
            if 'positive' in ds[zVarName].attrs and \
                    ds[zVarName].attrs['positive'] == 'down':
                attrs = ds[zVarName].attrs
                ds[zVarName] = -ds[zVarName]
                ds[zVarName].attrs = attrs
                ds[zVarName].attrs['positive'] = 'up'

            T, S, z = xarray.broadcast(ds[TVarName], ds[SVarName],
                                       ds[zVarName])

            ds['zBroadcast'] = z

            self.logger.info('  computing the dataset')
            ds.compute()
            self.logger.info('  writing {}...'.format(self.fileName))
            write_netcdf_with_fill(ds, self.fileName)
            for file_name in temp_files:
                self.logger.info('  Deleting temp file {}'.format(file_name))
                os.remove(file_name)
            self.logger.info('  Done!')

    def _get_file_name(self, obsDict, suffix=''):
        obsSection = '{}Observations'.format(self.componentName)
        climatologyDirectory = build_config_full_path(
            config=self.config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection=obsSection)

        make_directories(climatologyDirectory)

        fileName = '{}/{}_{}_{}{}.nc'.format(
            climatologyDirectory, 'TS_{}'.format(obsDict['suffix']),
            obsDict['gridName'], self.season, suffix)
        return fileName


class ComputeRegionTSSubtask(AnalysisTask):
    """
    Plots a T-S diagram for a given ocean region

    Attributes
    ----------
    regionGroup : str
        Name of the collection of region to plot

    regionName : str
        Name of the region to plot

    sectionName : str
        The section of the config file to get options from

    controlConfig : mpas_tools.config.MpasConfigParser
        The configuration options for the control run (if any)

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    mpasMasksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask MPAS files for each region to plot, used
        to get the mask file name

    obsDicts : dict of dicts
        Information on the observations to compare against

    season : str
        The season to compute the climatology for
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, regionName, controlConfig,
                 sectionName, fullSuffix, mpasClimatologyTask,
                 mpasMasksSubtask, obsDicts, season):

        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            Name of the collection of region to plot

        regionName : str
            Name of the region to plot

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)

        sectionName : str
            The config section with options for this regionGroup

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        mpasMasksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask MPAS files for each region to plot, used
            to get the mask file name

        obsDicts : dict of dicts
            Information on the observations to compare agains

        season : str
            The season to comput the climatogy for
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeRegionTSSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='compute{}_{}'.format(fullSuffix, season))

        self.run_after(mpasClimatologyTask)
        self.regionGroup = regionGroup
        self.regionName = regionName
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.mpasClimatologyTask = mpasClimatologyTask
        self.mpasMasksSubtask = mpasMasksSubtask
        self.obsDicts = obsDicts
        self.season = season
        self.prefix = fullSuffix

        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))

    def run_task(self):
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nComputing TS for {}...".format(self.regionName))

        zmin, zmax = self._write_mpas_t_s(self.config)

        for obsName in self.obsDicts:
            self._write_obs_t_s(self.obsDicts[obsName], zmin, zmax)

    def _write_mpas_t_s(self, config):

        climatologyName = 'TS_{}_{}'.format(self.prefix, self.season)
        outFileName = get_masked_mpas_climatology_file_name(
            config, self.season, self.componentName, climatologyName, op='avg')

        if os.path.exists(outFileName):
            ds = xarray.open_dataset(outFileName)
            zmin, zmax = ds.zbounds.values
            return zmin, zmax

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            self.logger.info('  Extracting T and S in the region...')

            sectionName = self.sectionName

            cellsChunk = 32768
            chunk = {'nCells': cellsChunk}

            try:
                restartFileName = self.runStreams.readpath('restart')[0]
            except ValueError:
                raise IOError('No MPAS-O restart file found: need at least one'
                              ' restart file to plot T-S diagrams')
            dsRestart = xarray.open_dataset(restartFileName)
            dsRestart = dsRestart.isel(Time=0)
            if 'landIceMask' in dsRestart:
                landIceMask = dsRestart.landIceMask
            else:
                landIceMask = None
            dsRestart = dsRestart.chunk(chunk)

            regionMaskFileName = self.mpasMasksSubtask.maskFileName

            dsRegionMask = xarray.open_dataset(regionMaskFileName)

            maskRegionNames = decode_strings(dsRegionMask.regionNames)
            regionIndex = maskRegionNames.index(self.regionName)

            dsMask = dsRegionMask.isel(nRegions=regionIndex)

            cellMask = dsMask.regionCellMasks == 1
            if landIceMask is not None:
                # only the region outside of ice-shelf cavities
                cellMask = numpy.logical_and(cellMask, landIceMask == 0)

            if config.has_option(sectionName, 'zmin'):
                zmin = config.getfloat(sectionName, 'zmin')
            else:
                if 'zminRegions' in dsMask:
                    zmin = dsMask.zminRegions.values
                else:
                    zmin = dsMask.zmin.values

            if config.has_option(sectionName, 'zmax'):
                zmax = config.getfloat(sectionName, 'zmax')
            else:
                if 'zmaxRegions' in dsMask:
                    zmax = dsMask.zmaxRegions.values
                else:
                    zmax = dsMask.zmax.values

            inFileName = get_unmasked_mpas_climatology_file_name(
                config, self.season, self.componentName, op='avg')

            ds = xarray.open_dataset(inFileName)

            variableList = ['timeMonthly_avg_activeTracers_temperature',
                            'timeMonthly_avg_activeTracers_salinity',
                            'timeMonthly_avg_layerThickness']
            ds = ds[variableList]

            ds['zMid'] = compute_zmid(dsRestart.bottomDepth,
                                      dsRestart.maxLevelCell-1,
                                      dsRestart.layerThickness)

            ds['volume'] = (dsRestart.areaCell *
                            ds['timeMonthly_avg_layerThickness'])

            ds.load()

            ds = ds.where(cellMask, drop=True)

            depthMask = numpy.logical_and(ds.zMid >= zmin,
                                          ds.zMid <= zmax)
            depthMask.compute()
            ds['depthMask'] = depthMask

            for var in variableList:
                ds[var] = ds[var].where(depthMask)

            T = ds['timeMonthly_avg_activeTracers_temperature'].values.ravel()
            mask = numpy.isfinite(T)
            T = T[mask]

            S = ds['timeMonthly_avg_activeTracers_salinity'].values.ravel()
            S = S[mask]

            zMid = ds['zMid'].values.ravel()[mask]

            volume = ds['volume'].values.ravel()[mask]

            dsOut = xarray.Dataset()
            dsOut['T'] = ('nPoints', T)
            dsOut['S'] = ('nPoints', S)
            dsOut['z'] = ('nPoints', zMid)
            dsOut['volume'] = ('nPoints', volume)
            dsOut['zbounds'] = ('nBounds', [zmin, zmax])
            write_netcdf_with_fill(dsOut, outFileName)

        return zmin, zmax

    def _write_obs_t_s(self, obsDict, zmin, zmax):
        obsSection = '{}Observations'.format(self.componentName)
        climatologyDirectory = build_config_full_path(
            config=self.config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection=obsSection)

        outFileName = '{}/TS_{}_{}_{}.nc'.format(
            climatologyDirectory, obsDict['suffix'], self.prefix, self.season)

        if os.path.exists(outFileName):
            return

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            chunk = {obsDict['latVar']: 400,
                     obsDict['lonVar']: 400}

            regionMaskFileName = obsDict['maskTask'].maskFileName

            dsRegionMask = \
                xarray.open_dataset(regionMaskFileName).stack(
                        nCells=(obsDict['latVar'], obsDict['lonVar']))
            dsRegionMask = dsRegionMask.reset_index('nCells').drop_vars(
                [obsDict['latVar'], obsDict['lonVar']])
            if 'nCells' in dsRegionMask.data_vars:
                dsRegionMask = dsRegionMask.drop_vars(['nCells'])

            maskRegionNames = decode_strings(dsRegionMask.regionNames)
            regionIndex = maskRegionNames.index(self.regionName)

            dsMask = dsRegionMask.isel(nRegions=regionIndex)

            if 'regionMasks' in dsMask:
                # this is the name used by the mask creation tool in mpas_tools
                maskVar = 'regionMasks'
            elif 'regionCellMasks' in dsMask:
                # this is the name used in the old mask creation tool in
                # mpas-analysis
                maskVar = 'regionCellMasks'
            else:
                raise ValueError(f'The file {regionMaskFileName} doesn\'t '
                                 f'contain a mask variable: regionMasks or '
                                 f'regionCellMasks')

            cellMask = dsMask[maskVar] == 1

            TVarName = obsDict['TVar']
            SVarName = obsDict['SVar']
            zVarName = obsDict['zVar']
            volVarName = obsDict['volVar']

            obsFileName = obsDict['climatologyTask'][self.season].fileName
            ds = xarray.open_dataset(obsFileName, chunks=chunk)
            ds = ds.stack(nCells=(obsDict['latVar'], obsDict['lonVar']))
            ds = ds.reset_index('nCells').drop_vars(
                [obsDict['latVar'], obsDict['lonVar']])
            if 'nCells' in ds.data_vars:
               ds = ds.drop_vars(['nCells'])

            ds = ds.where(cellMask, drop=True)

            cellsChunk = 32768
            chunk = {'nCells': cellsChunk}

            ds = ds.chunk(chunk)

            depthMask = numpy.logical_and(ds[zVarName] >= zmin,
                                          ds[zVarName] <= zmax)
            ds = ds.where(depthMask)
            ds.compute()

            T = ds[TVarName].values.ravel()
            mask = numpy.isfinite(T)
            T = T[mask]

            S = ds[SVarName].values.ravel()[mask]
            z = ds['zBroadcast'].values.ravel()[mask]

            volume = ds[volVarName].values.ravel()[mask]

            dsOut = xarray.Dataset()
            dsOut['T'] = ('nPoints', T)
            dsOut['S'] = ('nPoints', S)
            dsOut['z'] = ('nPoints', z)
            dsOut['volume'] = ('nPoints', volume)
            dsOut['zbounds'] = ('nBounds', [zmin, zmax])
            write_netcdf_with_fill(dsOut, outFileName)


class PlotRegionTSDiagramSubtask(AnalysisTask):
    """
    Plots a T-S diagram for a given ocean region

    Attributes
    ----------
    regionGroup : str
        Name of the collection of region to plot

    regionName : str
        Name of the region to plot

    sectionName : str
        The section of the config file to get options from

    controlConfig : mpas_tools.config.MpasConfigParser
        The configuration options for the control run (if any)

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    mpasMasksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask MPAS files for each region to plot, used
        to get the mask file name

    obsDicts : dict of dicts
        Information on the observations to compare against

    season : str
        The season to compute the climatology for
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, regionGroup, regionName, controlConfig,
                 sectionName, fullSuffix, mpasClimatologyTask,
                 mpasMasksSubtask, obsDicts, season):

        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        regionGroup : str
            Name of the collection of region to plot

        regionName : str
            Name of the region to plot

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)

        sectionName : str
            The config section with options for this regionGroup

        fullSuffix : str
            The regionGroup and regionName combined and modified to be
            appropriate as a task or file suffix

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        mpasMasksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask MPAS files for each region to plot, used
            to get the mask file name

        obsDicts : dict of dicts
            Information on the observations to compare agains

        season : str
            The season to comput the climatogy for
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotRegionTSDiagramSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName='plot{}_{}'.format(fullSuffix, season))

        self.run_after(mpasClimatologyTask)
        self.regionGroup = regionGroup
        self.regionName = regionName
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.mpasClimatologyTask = mpasClimatologyTask
        self.mpasMasksSubtask = mpasMasksSubtask
        self.obsDicts = obsDicts
        self.season = season
        self.prefix = fullSuffix

        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(PlotRegionTSDiagramSubtask, self).setup_and_check()

        self.xmlFileNames = ['{}/TS_diagram_{}_{}.xml'.format(
            self.plotsDirectory, self.prefix, self.season)]
        return

    def run_task(self):
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("\nPlotting TS diagram for {}"
                         "...".format(self.regionName))

        register_custom_colormaps()

        config = self.config
        sectionName = self.sectionName

        startYear = self.mpasClimatologyTask.startYear
        endYear = self.mpasClimatologyTask.endYear

        regionMaskFile = self.mpasMasksSubtask.geojsonFileName

        fcAll = read_feature_collection(regionMaskFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.regionName:
                fc.add_feature(feature)
                break

        self.logger.info('  Make plots...')

        groupLink = 'tsDiag' + self.regionGroup.replace(' ', '')

        if config.has_option(sectionName, 'titleFontSize'):
            titleFontSize = config.getint(sectionName, 'titleFontSize')
        else:
            titleFontSize = config.getint('plot', 'titleFontSize')

        if config.has_option(sectionName, 'axisFontSize'):
            axisFontSize = config.getint(sectionName, 'axisFontSize')
        else:
            axisFontSize = config.getint('plot', 'axisFontSize')

        if config.has_option(sectionName, 'defaultFontSize'):
            defaultFontSize = config.getint(sectionName, 'defaultFontSize')
        else:
            defaultFontSize = config.getint('plot', 'defaultFontSize')
        matplotlib.rc('font', size=defaultFontSize)

        axis_font = {'size': axisFontSize}
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}

        nSubplots = 1 + len(self.obsDicts)
        if self.controlConfig is not None:
            nSubplots += 1

        if nSubplots == 4:
            nCols = 2
            nRows = 2
        else:
            nCols = min(nSubplots, 3)
            nRows = (nSubplots-1)//3 + 1

        axisIndices = numpy.reshape(numpy.arange(nRows*nCols),
                                    (nRows, nCols))[::-1, :].ravel()

        width = 3 + 4.5*nCols
        height = 2 + 4*nRows

        # noinspection PyTypeChecker
        fig, axarray = plt.subplots(nrows=nRows, ncols=nCols, sharey=True,
                                    figsize=(width, height))

        if nSubplots == 1:
            axarray = numpy.array(axarray)

        if nRows == 1:
            axarray = axarray.reshape((nRows, nCols))

        T, S, zMid, volume, zmin, zmax = self._get_mpas_t_s(self.config)
        mainRunName = config.get('runs', 'mainRunName')
        plotFields = [{'S': S, 'T': T, 'z': zMid, 'vol': volume,
                       'title': mainRunName}]

        if self.controlConfig is not None:
            T, S, zMid, volume, _, _ = self._get_mpas_t_s(
                self.controlConfig)
            controlRunName = self.controlConfig.get('runs', 'mainRunName')
            plotFields.append({'S': S, 'T': T, 'z': zMid, 'vol': volume,
                               'title': 'Control: {}'.format(controlRunName)})

        for obsName in self.obsDicts:
            obsT, obsS, obsZ, obsVol = self._get_obs_t_s(
                self.obsDicts[obsName])
            plotFields.append({'S': obsS, 'T': obsT, 'z': obsZ, 'vol': obsVol,
                               'title': obsName})

        Tbins = config.getexpression(sectionName, 'Tbins', use_numpyfunc=True)
        Sbins = config.getexpression(sectionName, 'Sbins', use_numpyfunc=True)

        normType = config.get(sectionName, 'normType')

        PT, SP = numpy.meshgrid(Tbins, Sbins)
        SA = gsw.SA_from_SP(SP, p=0., lon=0., lat=-75.)
        CT = gsw.CT_from_t(SA, PT, p=0.)

        neutralDensity = sigma0(SA, CT)
        rhoInterval = config.getfloat(sectionName, 'rhoInterval')
        contours = numpy.arange(23., 29.+rhoInterval, rhoInterval)

        diagramType = config.get(sectionName, 'diagramType')
        if diagramType not in ['volumetric', 'scatter']:
            raise ValueError('Unexpected diagramType {}'.format(diagramType))

        lastPanel = None
        if config.has_option(sectionName, 'volMin'):
            volMinMpas = config.getfloat(sectionName, 'volMin')
        else:
            volMinMpas = None
        if config.has_option(sectionName, 'volMax'):
            volMaxMpas = config.getfloat(sectionName, 'volMax')
        else:
            volMaxMpas = None
        for index in range(len(axisIndices)):
            panelIndex = axisIndices[index]

            row = nRows-1 - index//nCols
            col = numpy.mod(index, nCols)

            if panelIndex >= nSubplots:
                plt.delaxes(axarray[row, col])
                continue

            plt.sca(axarray[row, col])
            T = plotFields[index]['T']
            S = plotFields[index]['S']
            z = plotFields[index]['z']
            volume = plotFields[index]['vol']
            title = plotFields[index]['title']
            title = limit_title(title, max_title_length=60)

            CS = plt.contour(SP, PT, neutralDensity, contours, linewidths=1.,
                             colors='k', zorder=2)
            plt.clabel(CS, fontsize=12, inline=1, fmt='%4.2f')

            if diagramType == 'volumetric':
                lastPanel, volMin, volMax = \
                    self._plot_volumetric_panel(T, S, volume)

                if volMinMpas is None:
                    volMinMpas = volMin
                if volMaxMpas is None:
                    volMaxMpas = volMax
                if normType == 'linear':
                    norm = colors.Normalize(vmin=0., vmax=volMaxMpas)
                elif normType == 'log':
                    if volMinMpas is None or volMaxMpas is None:
                        norm = None
                    else:
                        norm = colors.LogNorm(vmin=volMinMpas, vmax=volMaxMpas)
                else:
                    raise ValueError('Unsupported normType {}'.format(normType))
                if norm is not None:
                    lastPanel.set_norm(norm)
            else:
                lastPanel = self._plot_scatter_panel(T, S, z, zmin, zmax)

            CTFreezing = freezing.CT_freezing(Sbins, 0, 1)
            PTFreezing = gsw.t_from_CT(
                gsw.SA_from_SP(Sbins, p=0., lon=0., lat=-75.),
                CTFreezing, p=0.)
            plt.plot(Sbins, PTFreezing, linestyle='--', linewidth=1.,
                     color='k')

            plt.ylim([Tbins[0], Tbins[-1]])
            plt.xlim([Sbins[0], Sbins[-1]])

            plt.xlabel('Salinity (PSU)', **axis_font)
            if col == 0:
                plt.ylabel(r'Potential temperature ($^\circ$C)', **axis_font)
            plt.title(title, **title_font)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        fig.subplots_adjust(right=0.91)
        if nRows == 1:
            fig.subplots_adjust(top=0.85)
        else:
            fig.subplots_adjust(top=0.88)

        suptitle = 'T-S diagram for {} ({}, {:04d}-{:04d})\n' \
                   ' {} m < z < {} m'.format(self.regionName, self.season,
                                             startYear, endYear, zmin, zmax)
        fig.text(0.5, 0.9, suptitle, horizontalalignment='center',
                 **title_font)

        inset = add_inset(fig, fc, width=1.5, height=1.5)

        # add an empty plot covering the subplots to give common axis labels
        pos0 = axarray[0, 0].get_position()
        pos1 = axarray[-1, -1].get_position()
        pos_common = [pos0.x0, pos1.y0, pos1.x1-pos0.x0, pos0.y1-pos1.y0]
        common_ax = fig.add_axes(pos_common, zorder=-2)
        common_ax.spines['top'].set_color('none')
        common_ax.spines['bottom'].set_color('none')
        common_ax.spines['left'].set_color('none')
        common_ax.spines['right'].set_color('none')
        common_ax.tick_params(labelcolor='w', top=False, bottom=False,
                              left=False, right=False)

        common_ax.set_xlabel('Salinity (PSU)', **axis_font)
        common_ax.set_ylabel(r'Potential temperature ($^\circ$C)', **axis_font)

        # turn off labels for individual plots (just used for spacing)
        for index in range(len(axisIndices)):
            row = nRows-1 - index//nCols
            col = numpy.mod(index, nCols)
            ax = axarray[row, col]
            ax.set_xlabel('')
            ax.set_ylabel('')

        # move the color bar down a little ot avoid the inset
        pos0 = inset.get_position()
        pos1 = axarray[-1, -1].get_position()
        pad = 0.04
        top = pos0.y0 - pad
        height = top - pos1.y0
        cbar_ax = fig.add_axes([0.92, pos1.y0, 0.02, height])
        cbar = fig.colorbar(lastPanel, cax=cbar_ax)

        if diagramType == 'volumetric':
            cbar.ax.get_yaxis().labelpad = 15
            cbar.ax.set_ylabel(r'volume (m$^3$)', rotation=270, **axis_font)
        else:
            cbar.ax.set_ylabel('depth (m)', rotation=270, **axis_font)

        outFileName = '{}/TS_diagram_{}_{}.png'.format(
            self.plotsDirectory, self.prefix, self.season)
        savefig(outFileName, config, tight=False)

        caption = 'Regional mean of {}'.format(suptitle)
        write_image_xml(
            config=config,
            filePrefix='TS_diagram_{}_{}'.format(self.prefix, self.season),
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='T-S Diagrams',
            groupLink=groupLink,
            gallery=self.regionGroup,
            thumbnailDescription=self.regionName,
            imageDescription=caption,
            imageCaption=caption)

    def _get_mpas_t_s(self, config):
        climatologyName = 'TS_{}_{}'.format(self.prefix, self.season)
        inFileName = get_masked_mpas_climatology_file_name(
            config, self.season, self.componentName, climatologyName, op='avg')

        ds = xarray.open_dataset(inFileName)
        T = ds.T.values
        S = ds.S.values
        z = ds.z.values
        volume = ds.volume.values
        zmin, zmax = ds.zbounds.values

        return T, S, z, volume, zmin, zmax

    def _get_obs_t_s(self, obsDict):

        obsSection = '{}Observations'.format(self.componentName)
        climatologyDirectory = build_config_full_path(
            config=self.config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection=obsSection)

        inFileName = '{}/TS_{}_{}_{}.nc'.format(
            climatologyDirectory, obsDict['suffix'], self.prefix, self.season)

        ds = xarray.open_dataset(inFileName)
        T = ds.T.values
        S = ds.S.values
        z = ds.z.values
        volume = ds.volume.values

        return T, S, z, volume

    def _plot_volumetric_panel(self, T, S, volume):

        config = self.config
        sectionName = self.sectionName
        cmap = config.get(sectionName, 'colorMap')
        Tbins = config.getexpression(sectionName, 'Tbins',
                                     use_numpyfunc=True)
        Sbins = config.getexpression(sectionName, 'Sbins',
                                     use_numpyfunc=True)

        hist, _, _, panel = plt.hist2d(S, T, bins=[Sbins, Tbins],
                                       weights=volume, cmap=cmap, zorder=1,
                                       rasterized=True)

        poshist = hist[hist > 0.]
        if len(poshist) > 0:
            volMin = numpy.amin(poshist)
            volMax = numpy.amax(poshist)
        else:
            volMin = None
            volMax = None
        return panel, volMin, volMax

    def _plot_scatter_panel(self, T, S, z, zmin, zmax):

        config = self.config
        sectionName = self.sectionName
        cmap = config.get(sectionName, 'colorMap')

        indices = numpy.argsort(z)[::-1]

        panel = plt.scatter(S[indices], T[indices], c=z[indices],
                            s=5, vmin=zmin, vmax=zmax, cmap=cmap, zorder=1)

        return panel
