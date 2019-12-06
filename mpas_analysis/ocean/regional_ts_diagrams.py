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
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray
import numpy
import matplotlib.pyplot as plt
from matplotlib import colors
import gsw
try:
    from gsw import freezing
except ImportError:
    from gsw.gibbs import freezing
try:
    from gsw.density import sigma0
except ImportError:
    from gsw.gibbs.density_enthalpy_ct_exact import sigma0_CT_exact as sigma0

import dask
import multiprocessing
from multiprocessing.pool import ThreadPool

from geometric_features import FeatureCollection, read_feature_collection

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.plot import savefig, add_inset

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.io.utility import decode_strings, get_region_mask, \
    build_obs_path, build_config_full_path, make_directories

from mpas_analysis.shared.html import write_image_xml

from mpas_analysis.shared.regions import ComputeRegionMasksSubtask, \
    get_feature_list

from mpas_analysis.ocean.utility import compute_zmid

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.climatology import compute_climatology, \
    get_unmasked_mpas_climatology_file_name

from mpas_analysis.shared.plot.colormap import register_custom_colormaps


class RegionalTSDiagrams(AnalysisTask):  # {{{
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

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig :  ``MpasAnalysisConfigParser``, optional
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
            tags=['climatology', 'regions', 'publicObs'])

        self.run_after(mpasClimatologyTask)
        self.mpasClimatologyTask = mpasClimatologyTask

        regionGroups = config.getExpression(self.taskName, 'regionGroups')

        self.seasons = config.getExpression(self.taskName, 'seasons')

        parallelTaskCount = config.getWithDefault('execute',
                                                  'parallelTaskCount',
                                                  default=1)

        obsDicts = {'SOSE': {
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
            'lonVar': 'lon',
            'latVar': 'lat',
            'TVar': 'theta',
            'SVar': 'salinity',
            'volVar': 'volume',
            'zVar': 'z',
            'tVar': 'Time'}, 'WOA18': {
            'suffix': 'WOA18',
            'gridName': 'Global_0.25x0.25degree',
            'gridFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'TFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'SFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'volFileName': None,
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
            obsList = config.getExpression(sectionName, 'obs')
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
            fileSuffix = sectionSuffix[0].lower() + sectionSuffix[1:]

            regionMaskFile = config.getExpression(sectionName,
                                                  'regionMask')

            regionMaskFile = get_region_mask(config, regionMaskFile)

            regionNames = config.getExpression(sectionName, 'regionNames')

            if 'all' in regionNames and os.path.exists(regionMaskFile):
                regionNames = get_feature_list(regionMaskFile)

            subtaskName = 'computeMpas{}Masks'.format(sectionSuffix)
            mpasMasksSubtask = ComputeRegionMasksSubtask(
                self, regionMaskFile, outFileSuffix=fileSuffix,
                featureList=regionNames, subtaskName=subtaskName,
                subprocessCount=parallelTaskCount)

            self.add_subtask(mpasMasksSubtask)

            obsList = config.getExpression(sectionName, 'obs')
            groupObsDicts = {}

            for obsName in obsList:
                localObsDict = dict(obsDicts[obsName])
                suffix = localObsDict['suffix']
                subtaskName = 'compute{}{}Masks'.format(suffix, sectionSuffix)
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=localObsDict['gridFileName'])
                obsMasksSubtask = ComputeRegionMasksSubtask(
                    self, regionMaskFile, outFileSuffix=fileSuffix,
                    featureList=regionNames, subtaskName=subtaskName,
                    subprocessCount=parallelTaskCount,
                    obsFileName=obsFileName, lonVar=localObsDict['lonVar'],
                    latVar=localObsDict['latVar'],
                    meshName=localObsDict['gridName'])

                obsDicts[obsName]['maskTask'] = obsMasksSubtask

                self.add_subtask(obsMasksSubtask)
                localObsDict['maskTask'] = obsMasksSubtask
                groupObsDicts[obsName] = localObsDict

            for regionName in regionNames:

                fullSuffix = sectionSuffix + '_' + \
                    regionName[0].lower() + \
                    regionName[1:].replace(' ', '')

                for season in self.seasons:
                    plotRegionSubtask = PlotRegionTSDiagramSubtask(
                        self, regionGroup, regionName, controlConfig,
                        sectionName, fullSuffix, mpasClimatologyTask,
                        mpasMasksSubtask, groupObsDicts, season)
                    for obsName in obsList:
                        task = \
                            obsDicts[obsName]['climatologyTask'][season]
                        plotRegionSubtask.run_after(task)
                        task = obsDicts[obsName]['maskTask']
                        plotRegionSubtask.run_after(task)
                    self.add_subtask(plotRegionSubtask)

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
        super(RegionalTSDiagrams, self).setup_and_check()

        # don't add the variables and seasons to mpasClimatologyTask until
        # we're sure this subtask is supposed to run
        variableList = ['timeMonthly_avg_activeTracers_temperature',
                        'timeMonthly_avg_activeTracers_salinity',
                        'timeMonthly_avg_layerThickness']
        self.mpasClimatologyTask.add_variables(variableList=variableList,
                                               seasons=self.seasons)
    # }}}


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

    def __init__(self, parentTask, obsName, obsDict, season):  # {{{
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

        # }}}

    def run_task(self):  # {{{
        """
        Plots time-series output of properties in an ocean region.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if os.path.exists(self.fileName):
            return

        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            self.logger.info("\n computing T S climatogy for {}...".format(
                self.obsName))

            config = self.config
            obsDict = self.obsDict

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
            ds = xarray.open_dataset(obsFileName, chunks=chunk)
            if obsDict['SFileName'] != obsDict['TFileName']:
                obsFileName = build_obs_path(
                    config, component=self.componentName,
                    relativePath=obsDict['SFileName'])
                dsS = xarray.open_dataset(obsFileName, chunks=chunk)
                ds[SVarName] = dsS[SVarName]

            if obsDict['volFileName'] is None:
                # compute volume from lat, lon, depth bounds
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
                dsVol = xarray.open_dataset(obsFileName)
                ds[volVarName] = dsVol[volVarName]

            ds.compute()

            chunk = {obsDict['latVar']: 400,
                     obsDict['lonVar']: 400}

            ds = ds.chunk(chunk)

            if obsDict['tVar'] in ds.dims:
                if obsDict['tVar'] != 'Time':
                    if obsDict['tVar'] == 'month':
                        ds = ds.rename({obsDict['tVar']: 'Time'})
                        ds.coords['month'] = ds['Time']
                    else:
                        ds = ds.rename({obsDict['tVar']: 'Time'})
                if 'year' not in ds:
                    ds.coords['year'] = numpy.ones(ds.sizes['Time'], int)

                monthValues = constants.monthDictionary[self.season]
                ds = compute_climatology(ds, monthValues, maskVaries=True)

            if 'positive' in ds[zVarName].attrs and \
                    ds[zVarName].attrs['positive'] == 'down':
                attrs = ds[zVarName].attrs
                ds[zVarName] = -ds[zVarName]
                ds[zVarName].attrs = attrs
                ds[zVarName].attrs['positive'] = 'up'

            T, S, z = xarray.broadcast(ds[TVarName], ds[SVarName],
                                       ds[zVarName])

            ds['zBroadcast'] = z

            write_netcdf(ds, self.fileName)  # }}}

    def _get_file_name(self, obsDict):
        obsSection = '{}Observations'.format(self.componentName)
        climatologyDirectory = build_config_full_path(
            config=self.config, section='output',
            relativePathOption='climatologySubdirectory',
            relativePathSection=obsSection)

        make_directories(climatologyDirectory)

        fileName = '{}/{}_{}_{}.nc'.format(
            climatologyDirectory, 'TS_{}'.format(obsDict['suffix']),
            obsDict['gridName'], self.season)
        return fileName

    # }}}


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

    controlConfig : ``MpasAnalysisConfigParser``
        The configuration options for the control run (if any)

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

    def __init__(self, parentTask, regionGroup, regionName, controlConfig,
                 sectionName, fullSuffix, mpasClimatologyTask,
                 mpasMasksSubtask, obsDicts, season):
        # {{{
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

        controlConfig :  ``MpasAnalysisConfigParser``, optional
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

        self.regionGroup = regionGroup
        self.regionName = regionName
        self.sectionName = sectionName
        self.controlConfig = controlConfig
        self.mpasClimatologyTask = mpasClimatologyTask
        self.mpasMasksSubtask = mpasMasksSubtask
        self.obsDicts = obsDicts
        self.season = season
        self.prefix = fullSuffix[0].lower() + fullSuffix[1:]

        parallelTaskCount = self.config.getint('execute', 'parallelTaskCount')
        self.subprocessCount = min(parallelTaskCount,
                                   self.config.getint(self.taskName,
                                                      'subprocessCount'))
        self.daskThreads = min(
            multiprocessing.cpu_count(),
            self.config.getint(self.taskName, 'daskThreads'))
        # }}}

    def setup_and_check(self):  # {{{
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

        self.xmlFileNames = ['{}/{}_{}.xml'.format(
            self.plotsDirectory, self.prefix, self.season)]
        return  # }}}

    def run_task(self):  # {{{
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

        regionMaskFile = config.getExpression(sectionName, 'regionMask')

        regionMaskFile = get_region_mask(config, regionMaskFile)

        fcAll = read_feature_collection(regionMaskFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.regionName:
                fc.add_feature(feature)
                break

        self.logger.info('  Make plots...')

        groupLink = 'tsDiag' + self.regionGroup[0].lower() + \
            self.regionGroup[1:].replace(' ', '')

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

        titleFontSize = config.get('plot', 'titleFontSize')
        axis_font = {'size': config.get('plot', 'axisFontSize')}
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}

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
                self.obsDicts[obsName], zmin, zmax)
            plotFields.append({'S': obsS, 'T': obsT, 'z': obsZ, 'vol': obsVol,
                               'title': obsName})

        Tbins = config.getExpression(sectionName, 'Tbins', usenumpyfunc=True)
        Sbins = config.getExpression(sectionName, 'Sbins', usenumpyfunc=True)

        normType = config.get(sectionName, 'normType')

        PT, SP = numpy.meshgrid(Tbins, Sbins)
        SA = gsw.SA_from_SP(SP, p=0., lon=0., lat=-75.)
        CT = gsw.CT_from_t(SA, PT, p=0.)

        neutralDensity = sigma0(SA, CT)
        rhoInterval = config.getfloat(sectionName, 'rhoInterval')
        contours = numpy.arange(24., 29.+rhoInterval, rhoInterval)

        diagramType = config.get(sectionName, 'diagramType')
        if diagramType not in ['volumetric', 'scatter']:
            raise ValueError('Unexpected diagramType {}'.format(diagramType))

        lastPanel = None
        volMinMpas = None
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

            CS = plt.contour(SP, PT, neutralDensity, contours, linewidths=1.,
                             colors='k', zorder=2)
            plt.clabel(CS, fontsize=12, inline=1, fmt='%4.2f')

            if diagramType == 'volumetric':
                lastPanel, volMin, volMax = \
                    self._plot_volumetric_panel(T, S, volume)

                if index == 0:
                    volMinMpas = volMin
                    volMaxMpas = volMax
                if normType == 'linear':
                    norm = colors.Normalize(vmin=0., vmax=volMaxMpas)
                elif normType == 'log':
                    norm = colors.LogNorm(vmin=volMinMpas, vmax=volMaxMpas)
                else:
                    raise ValueError('Unsupported normType {}'.format(normType))
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
            plt.title(title)

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
            cbar.ax.set_ylabel(r'volume (m$^3$)', rotation=270)
        else:
            cbar.ax.set_ylabel('depth (m)', rotation=270)

        outFileName = '{}/{}_{}.png'.format(self.plotsDirectory, self.prefix,
                                            self.season)
        savefig(outFileName, tight=False)

        caption = 'Regional mean of {}'.format(suptitle)
        write_image_xml(
            config=config,
            filePrefix='{}_{}'.format(self.prefix, self.season),
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='T-S Diagrams',
            groupLink=groupLink,
            gallery=self.regionGroup,
            thumbnailDescription=self.regionName,
            imageDescription=caption,
            imageCaption=caption)

        # }}}

    def _get_mpas_t_s(self, config):  # {{{
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
            dsRestart = dsRestart.isel(Time=0).chunk(chunk)

            regionMaskFileName = self.mpasMasksSubtask.maskFileName

            dsRegionMask = xarray.open_dataset(regionMaskFileName)

            maskRegionNames = decode_strings(dsRegionMask.regionNames)
            regionIndex = maskRegionNames.index(self.regionName)

            dsMask = dsRegionMask.isel(nRegions=regionIndex).chunk(chunk)

            cellMask = dsMask.regionCellMasks == 1
            if 'landIceMask' in dsRestart:
                # only the region outside of ice-shelf cavities
                cellMask = numpy.logical_and(cellMask,
                                             dsRestart.landIceMask == 0)

            if config.has_option(sectionName, 'zmin'):
                zmin = config.getfloat(sectionName, 'zmin')
            else:
                zmin = dsMask.zmin.values

            if config.has_option(sectionName, 'zmax'):
                zmax = config.getfloat(sectionName, 'zmax')
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
                                      dsRestart.maxLevelCell,
                                      dsRestart.layerThickness)

            ds['volume'] = (dsRestart.areaCell *
                            ds['timeMonthly_avg_layerThickness'])

            ds = ds.where(cellMask, drop=True)

            self.logger.info("Don't worry about the following dask "
                             "warnings.")
            depthMask = numpy.logical_and(ds.zMid >= zmin,
                                          ds.zMid <= zmax)
            depthMask.compute()
            self.logger.info("Dask warnings should be done.")
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

        return T, S, zMid, volume, zmin, zmax  # }}}

    def _get_obs_t_s(self, obsDict, zmin, zmax):  # {{{
        with dask.config.set(schedular='threads',
                             pool=ThreadPool(self.daskThreads)):

            chunk = {obsDict['latVar']: 400,
                     obsDict['lonVar']: 400}

            regionMaskFileName = obsDict['maskTask'].maskFileName

            dsRegionMask = \
                xarray.open_dataset(regionMaskFileName).chunk(chunk).stack(
                        nCells=(obsDict['latVar'], obsDict['lonVar']))
            dsRegionMask = dsRegionMask.reset_index('nCells').drop_vars(
                [obsDict['latVar'], obsDict['lonVar']])

            maskRegionNames = decode_strings(dsRegionMask.regionNames)
            regionIndex = maskRegionNames.index(self.regionName)

            dsMask = dsRegionMask.isel(nRegions=regionIndex)

            cellMask = dsMask.regionCellMasks == 1

            TVarName = obsDict['TVar']
            SVarName = obsDict['SVar']
            zVarName = obsDict['zVar']
            volVarName = obsDict['volVar']

            obsFileName = obsDict['climatologyTask'][self.season].fileName
            ds = xarray.open_dataset(obsFileName, chunks=chunk)
            ds = ds.stack(nCells=(obsDict['latVar'], obsDict['lonVar']))
            ds = ds.reset_index('nCells').drop_vars(
                [obsDict['latVar'], obsDict['lonVar']])

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

        return T, S, z, volume  # }}}

    def _plot_volumetric_panel(self, T, S, volume):  # {{{

        config = self.config
        sectionName = self.sectionName
        cmap = config.get(sectionName, 'colorMap')
        Tbins = config.getExpression(sectionName, 'Tbins',
                                     usenumpyfunc=True)
        Sbins = config.getExpression(sectionName, 'Sbins',
                                     usenumpyfunc=True)

        hist, _, _, panel = plt.hist2d(S, T, bins=[Sbins, Tbins],
                                       weights=volume, cmap=cmap, zorder=1)

        poshist = hist[hist > 0.]
        if len(poshist) > 0:
            volMin = numpy.amin(poshist)
            volMax = numpy.amax(poshist)
        else:
            volMin = None
            volMax = None
        return panel, volMin, volMax  # }}}

    def _plot_scatter_panel(self, T, S, z, zmin, zmax):  # {{{

        config = self.config
        sectionName = self.sectionName
        cmap = config.get(sectionName, 'colorMap')

        indices = numpy.argsort(z)[::-1]

        panel = plt.scatter(S[indices], T[indices], c=z[indices],
                            s=5, vmin=zmin, vmax=zmax, cmap=cmap, zorder=1)

        return panel  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
