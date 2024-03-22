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
import csv

import numpy as np
import xarray as xr
import dask
from multiprocessing.pool import ThreadPool
from pyremap import ProjectionGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path, decode_strings, \
    build_config_full_path
from mpas_analysis.shared.io import write_netcdf_with_fill

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask
from mpas_analysis.shared.projection import get_pyproj_projection
from mpas_analysis.shared.climatology.climatology import \
    get_masked_mpas_climatology_file_name

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.constants import constants


class ClimatologyMapAntarcticMelt(AnalysisTask):
    """
    An analysis task for comparison of Antarctic melt rates against
    observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, regionMasksTask,
                 controlConfig):
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

        controlConfig : mpas_tools.config.MpasConfigParser
            Configuration options for a control run
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'meltRate'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapAntarcticMelt, self).__init__(
            config=config, taskName='climatologyMapAntarcticMelt',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName,
                  'landIceCavities', 'antarctic'])

        sectionName = self.taskName

        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             f'valid list of seasons')

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        makeTables = config.getboolean(sectionName, 'makeTables')

        if makeTables:
            for season in seasons:
                tableSubtask = AntarcticMeltTableSubtask(
                    parentTask=self, mpasClimatologyTask=mpasClimatologyTask,
                    controlConfig=controlConfig,
                    regionMasksTask=regionMasksTask, season=season)
                self.add_subtask(tableSubtask)

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             f'valid list of comparison grids')

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasAntarcticMeltClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=None,
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:

            refTitleLabel = \
                'Observations (Paolo et al. 2023)'

            observationsDirectory = build_obs_path(
                config, 'ocean', 'meltSubdirectory')

            comparison_res = config.getfloat(
                'climatology', 'comparisonAntarcticStereoResolution')

            # the maximum available resolution that is not coarser than the
            # comparison
            avail_res = np.array([10., 4., 1.])
            valid = avail_res >= comparison_res
            if np.count_nonzero(valid) == 0:
                res = np.amin(avail_res)
            else:
                res = np.amax(avail_res[valid])

            obsFileName = \
                f'{observationsDirectory}/Paolo/Paolo_2023_' \
                f'iceshelf_melt_rates_1992-2017_v1.0_6000x6000km_{res:g}km_' \
                f'Antarctic_stereo.20240220.nc'

            remapObservationsSubtask = RemapObservedAntarcticMeltClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix='meltRate',
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'
            diffTitleLabel = 'Main - Control'

        totalFluxVar = 'timeMonthly_avg_landIceFreshwaterFluxTotal'
        landIceFluxVar = 'timeMonthly_avg_landIceFreshwaterFlux'
        frazilFluxVar = 'timeMonthly_avg_frazilIceFreshwaterFlux'

        mpasFieldName = totalFluxVar

        if controlConfig is None:
            refFieldName = 'meltRate'
        else:
            refFieldName = mpasFieldName

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtaskName = f'plot_total_melt_{season}_{comparisonGridName}'
                subtask = PlotAntarcticMeltSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig,
                    subtaskName=subtaskName)

                subtask.set_plot_info(
                    outFileLabel='antMeltTotal',
                    fieldNameInTitle='Total Melt Flux',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm a$^{-1}$ freshwater equiv.',
                    imageCaption='Antarctic Total Melt Flux',
                    galleryGroup='Melt Rate',
                    groupSubtitle=None,
                    groupLink='antarctic_melt',
                    galleryName='Total Melt Flux')

                self.add_subtask(subtask)

        mpasFieldName = landIceFluxVar

        if controlConfig is None:
            refFieldName = 'meltRate'
        else:
            refFieldName = mpasFieldName

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtaskName = \
                    f'plot_interface_melt_{season}_{comparisonGridName}'
                subtask = PlotAntarcticMeltSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig,
                    subtaskName=subtaskName)

                # In PlotAntarcticMeltSubtask, we will remove the obs from
                # these plots if totalFluxVar is present so we only compare one
                # field with obs

                subtask.set_plot_info(
                    outFileLabel='antMeltInterface',
                    fieldNameInTitle='Melt Rate at Interface',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm a$^{-1}$ freshwater equiv.',
                    imageCaption='Antarctic Melt Rate at Interface',
                    galleryGroup='Melt Rate',
                    groupSubtitle=None,
                    groupLink='antarctic_melt_int',
                    galleryName='Melt Rate at the Ice-ocean Interface')

                self.add_subtask(subtask)

        mpasFieldName = frazilFluxVar

        if controlConfig is None:
            refTitleLabel = None
            refFieldName = None
            diffTitleLabel = None

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'
            refFieldName = mpasFieldName
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtaskName = \
                    f'plot_interface_frazil_{season}_{comparisonGridName}'
                subtask = PlotAntarcticMeltSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    controlConfig=controlConfig, subtaskName=subtaskName)

                subtask.set_plot_info(
                    outFileLabel='antFrazil',
                    fieldNameInTitle='Frazil Accretion Rate, neg. upward',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm a$^{-1}$ freshwater equiv.',
                    imageCaption='Antarctic Accretion Rate',
                    galleryGroup='Melt Rate',
                    groupSubtitle=None,
                    groupLink='antarctic_frazil_flux',
                    galleryName='Frazil Accretion Rate')

                self.add_subtask(subtask)

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class
        # (AnalysisTask), which will perform some common setup
        super(ClimatologyMapAntarcticMelt, self).setup_and_check()

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['data', 'standalone', 'coupled']:
            raise ValueError('*** climatologyMapMeltAntarctic requires '
                             'config_land_ice_flux_mode \n'
                             '    to be data, standalone or coupled. '
                             '    Otherwise, no melt rates are available \n'
                             '    for plotting.')


class RemapMpasAntarcticMeltClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for remapping climatologies of Antarctic melt rates and adding

    Attributes
    ----------
    landIceMask : xarray.DataArray
        A mask indicating where there is land ice on the ocean grid (thus,
        where melt rates are valid)

    renameDict : dict
        A dictionary use to rename variables in the climatology
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def setup_and_check(self):
        """
        Figure out which variable(s) to remap
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        totalFluxVar = 'timeMonthly_avg_landIceFreshwaterFluxTotal'
        landIceFluxVar = 'timeMonthly_avg_landIceFreshwaterFlux'
        frazilFluxVar = 'timeMonthly_avg_frazilIceFreshwaterFlux'

        if totalFluxVar in self.mpasClimatologyTask.allVariables:
            # include the total and constituent fluxes
            self.variableList = [totalFluxVar, landIceFluxVar, frazilFluxVar]
        else:
            # we only have the old name without the frazil accretion rate
            self.variableList = [landIceFluxVar]

        super().setup_and_check()

    def run_task(self):
        """
        Compute climatologies of melt rates from E3SM/MPAS output

        This function has been overridden to load ``landIceMask`` from a
        restart file for later use in masking the melt rate.  It then simply
        calls the run function from
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, load the land-ice mask from the restart file
        dsLandIceMask = xr.open_dataset(self.restartFileName)
        dsLandIceMask = dsLandIceMask[['landIceMask']]
        dsLandIceMask = dsLandIceMask.isel(Time=0)
        self.landIceMask = dsLandIceMask.landIceMask > 0.

        # then, call run from the base class (RemapMpasClimatologySubtask),
        # which will perform the main function of the task
        super(RemapMpasAntarcticMeltClimatology, self).run_task()

    def customize_masked_climatology(self, climatology, season):
        """
        Mask the melt rates using ``landIceMask`` and rescale it to m/yr

        Parameters
        ----------
        climatology : ``xarray.Dataset`` object
            the climatology data set

        season : str
            The name of the season to be masked

        Returns
        -------
        climatology : ``xarray.Dataset`` object
            the modified climatology data set
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        for fieldName in self.variableList:

            # scale the field to m/yr from kg/m^2/s and mask out non-land-ice
            # areas
            climatology[fieldName] = \
                constants.sec_per_year / constants.rho_fw * \
                climatology[fieldName].where(self.landIceMask)
            climatology[fieldName].attrs['units'] = 'm yr^-1'

        return climatology


class RemapObservedAntarcticMeltClimatology(RemapObservedClimatologySubtask):

    """
    A subtask for reading and remapping Antarctic melt-rate observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # create a descriptor of the observation grid using the x/y polar
        # stereographic coordinates
        projection = get_pyproj_projection(comparison_grid_name='antarctic')
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=fileName, xVarName='x', yVarName='y')

        # update the mesh name to match the format used elsewhere in
        # MPAS-Analysis
        x = obsDescriptor.x
        y = obsDescriptor.y
        width = 1e-3 * (x[-1] - x[0])
        height = 1e-3 * (y[-1] - y[0])
        res = 1e-3 * (x[1] - x[0])
        obsDescriptor.meshName = f'{width}x{height}km_{res}km_Antarctic_stereo'

        return obsDescriptor

    def build_observational_dataset(self, fileName):
        """
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
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # Load MLD observational data
        dsObs = xr.open_dataset(fileName)

        return dsObs


class AntarcticMeltTableSubtask(AnalysisTask):
    def __init__(self, parentTask, mpasClimatologyTask, controlConfig,
                 regionMasksTask, season, subtaskName=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  ``ClimatologyMapAntarcticMelt``
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlConfig : mpas_tools.config.MpasConfigParser
            Configuration options for a control run (if any)

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        season : str
            One of the seasons in ``constants.monthDictionary``

        subtaskName : str, optional
            The name of the subtask
        """
        # Authors
        # -------
        # Xylar Asay-Davis
        tags = ['climatology', 'table']

        if subtaskName is None:
            subtaskName = f'table{season}'

        # call the constructor from the base class (AnalysisTask)
        super(AntarcticMeltTableSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=tags)

        self.season = season
        self.mpasClimatologyTask = mpasClimatologyTask
        self.controlConfig = controlConfig

        self.masksSubtask = regionMasksTask.add_mask_subtask(
            regionGroup='Ice Shelves')
        self.iceShelfMasksFile = self.masksSubtask.geojsonFileName

        self.run_after(self.masksSubtask)
        self.run_after(mpasClimatologyTask)

    def run_task(self):
        """
        Computes and plots table of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.logger.info("Computing Antarctic melt rate table...")
        config = self.config

        sectionName = self.taskName
        iceShelvesInTable = config.getexpression(sectionName,
                                                 'iceShelvesInTable')
        if len(iceShelvesInTable) == 0:
            return

        iceShelvesInTable = self.masksSubtask.expand_region_names(
            iceShelvesInTable)

        meltRateFileName = get_masked_mpas_climatology_file_name(
            config, self.season, self.componentName,
            climatologyName='antarcticMeltTable')

        if not os.path.exists(meltRateFileName):
            with dask.config.set(schedular='threads',
                                 pool=ThreadPool(1)):

                # Load data:
                inFileName = \
                    self.mpasClimatologyTask.get_file_name(self.season)
                mpasFieldName = 'timeMonthly_avg_landIceFreshwaterFlux'
                dsIn = xr.open_dataset(inFileName)
                freshwaterFlux = dsIn[mpasFieldName]
                if 'Time' in freshwaterFlux.dims:
                    freshwaterFlux = freshwaterFlux.isel(Time=0)

                regionMaskFileName = self.masksSubtask.maskFileName

                dsRegionMask = xr.open_dataset(regionMaskFileName)

                # figure out the indices of the regions to plot
                regionNames = decode_strings(dsRegionMask.regionNames)

                regionIndices = []
                for iceShelf in iceShelvesInTable:
                    for index, regionName in enumerate(regionNames):
                        if iceShelf == regionName:
                            regionIndices.append(index)
                            break

                # select only those regions we want to plot
                dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
                cellMasks = \
                    dsRegionMask.regionCellMasks.chunk({'nRegions': 10})

                restartFileName = \
                    self.runStreams.readpath('restart')[0]

                dsRestart = xr.open_dataset(restartFileName)
                landIceFraction = dsRestart.landIceFraction.isel(Time=0)
                areaCell = dsRestart.areaCell

                # convert from kg/s to kg/yr
                totalMeltFlux = constants.sec_per_year * \
                    (cellMasks * areaCell * freshwaterFlux).sum(dim='nCells')
                totalMeltFlux.compute()

                totalArea = \
                    (landIceFraction * cellMasks * areaCell).sum(dim='nCells')

                # from kg/m^2/yr to m/yr
                meltRates = ((1. / constants.rho_fw) *
                             (totalMeltFlux / totalArea))
                meltRates.compute()

                # convert from kg/yr to GT/yr
                totalMeltFlux /= constants.kg_per_GT

                ds = xr.Dataset()
                ds['totalMeltFlux'] = totalMeltFlux
                ds.totalMeltFlux.attrs['units'] = 'GT a$^{-1}$'
                ds.totalMeltFlux.attrs['description'] = \
                    'Total melt flux summed over each ice shelf or region'
                ds['meltRates'] = meltRates
                ds.meltRates.attrs['units'] = 'm a$^{-1}$'
                ds.meltRates.attrs['description'] = \
                    'Melt rate averaged over each ice shelf or region'

                ds['area'] = 1e-6*totalArea
                ds.meltRates.attrs['units'] = 'km$^2$'
                ds.meltRates.attrs['description'] = \
                    'Region or ice shelf area'

                ds['regionNames'] = dsRegionMask.regionNames

                write_netcdf_with_fill(ds, meltRateFileName)
        else:
            ds = xr.open_dataset(meltRateFileName)

        mainRunName = config.get('runs', 'mainRunName')
        fieldNames = ['Region', 'Area', mainRunName]

        controlConfig = self.controlConfig
        if controlConfig is not None:
            controlFileName = get_masked_mpas_climatology_file_name(
                controlConfig, self.season, self.componentName,
                climatologyName='antarcticMeltTable')
            dsControl = xr.open_dataset(controlFileName)
            controlRunName = controlConfig.get('runs', 'mainRunName')
            fieldNames.append(controlRunName)
        else:
            dsControl = None
            controlRunName = None

        regionNames = decode_strings(ds.regionNames)

        tableBase = build_config_full_path(config, 'output',
                                           'tablesSubdirectory')
        outDirectory = f'{tableBase}/antarcticMelt/'

        try:
            os.makedirs(outDirectory)
        except OSError:
            pass

        tableFileName = \
            f'{outDirectory}/antarcticMeltRateTable_{self.season}.csv'

        with open(tableFileName, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)

            writer.writeheader()
            for index, regionName in enumerate(regionNames):
                row = {'Region': regionName,
                       'Area': f'{ds.area[index].values}',
                       mainRunName: f'{ds.meltRates[index].values}'}
                if dsControl is not None:
                    row[controlRunName] = \
                        f'{dsControl.meltRates[index].values}'
                writer.writerow(row)

        tableFileName = \
            f'{outDirectory}/antarcticMeltFluxTable_{self.season}.csv'

        with open(tableFileName, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)

            writer.writeheader()
            for index, regionName in enumerate(regionNames):
                row = {'Region': regionName,
                       'Area': f'{ds.area[index].values}',
                       mainRunName: f'{ds.totalMeltFlux[index].values}'}
                if dsControl is not None:
                    row[controlRunName] = \
                        f'{dsControl.totalMeltFlux[index].values}'
                writer.writerow(row)


class PlotAntarcticMeltSubtask(PlotClimatologyMapSubtask):
    """
    A subtask for plotting antarctic melt fields if available

    Attributes
    ----------
    doPlot : bool
        Whether the required variable from the climatology is available so that
        a plot should be generated
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.
        """
        allVariables = \
            self.remapMpasClimatologySubtask.mpasClimatologyTask.allVariables

        totalFluxVar = 'timeMonthly_avg_landIceFreshwaterFluxTotal'
        landIceFluxVar = 'timeMonthly_avg_landIceFreshwaterFlux'
        plotAll = (totalFluxVar in allVariables)

        if self.mpasFieldName == landIceFluxVar and plotAll and \
                self.controlConfig is None:
            # need to remove obs because we only wnat to plot them vs the
            # total flux
            self.remapObsClimatologySubtask = None
            self.refTitleLabel = None
            self.refFieldName = None
            self.diffTitleLabel = None

        self.doPlot = (self.mpasFieldName == landIceFluxVar or plotAll)

        if self.doPlot:
            super().setup_and_check()
        else:
            # still need to call the base class's method
            AnalysisTask.setup_and_check(self=self)

    def run_task(self):
        """
        Plot the variable if available
        """
        if self.doPlot:
            super().run_task()
