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

        controlconfig : mpas_tools.config.MpasConfigParser
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

        mpasFieldName = 'timeMonthly_avg_landIceFreshwaterFlux'
        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

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
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapMpasAntarcticMeltClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:

            refTitleLabel = \
                'Observations (Rignot et al, 2013)'

            observationsDirectory = build_obs_path(
                config, 'ocean', 'meltSubdirectory')

            obsFileName = \
                '{}/Rignot_2013_melt_rates_6000.0x6000.0km_10.0km_' \
                'Antarctic_stereo.nc'.format(observationsDirectory)
            refFieldName = 'meltRate'
            outFileLabel = 'meltRignot'
            galleryName = 'Observations: Rignot et al. (2013)'

            remapObservationsSubtask = RemapObservedAntarcticMeltClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'melt'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='Melt Rate',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'm a$^{-1}$',
                    imageCaption='Antarctic Melt Rate',
                    galleryGroup='Melt Rate',
                    groupSubtitle=None,
                    groupLink='antarctic_melt',
                    galleryName=galleryName)

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
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** climatologyMapMeltAntarctic requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')


class RemapMpasAntarcticMeltClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for remapping climatologies of Antarctic melt rates and adding

    Attributes
    ----------
    landIceMask : xarray.DataArray
        A mask indicating where there is land ice on the ocean grid (thus,
        where melt rates are valid)
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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

        fieldName = self.variableList[0]

        # scale the field to m/yr from kg/m^2/s and mask out non-land-ice areas
        climatology[fieldName] = \
            constants.sec_per_year / constants.rho_fw * \
            climatology[fieldName].where(self.landIceMask)

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

        controlconfig : mpas_tools.config.MpasConfigParser
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
            subtaskName = 'table{}'.format(season)

        # call the constructor from the base class (AnalysisTask)
        super(AntarcticMeltTableSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            subtaskName=subtaskName,
            componentName=parentTask.componentName,
            tags=tags)

        config = parentTask.config
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
                inFileName = self.mpasClimatologyTask.get_file_name(self.season)
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
                cellMasks = dsRegionMask.regionCellMasks.chunk({'nRegions': 10})

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

        outDirectory = '{}/antarcticMelt/'.format(
            build_config_full_path(config, 'output', 'tablesSubdirectory'))

        try:
            os.makedirs(outDirectory)
        except OSError:
            pass

        tableFileName = '{}/antarcticMeltRateTable_{}.csv'.format(outDirectory,
                                                                  self.season)

        with open(tableFileName, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)

            writer.writeheader()
            for index, regionName in enumerate(regionNames):
                row = {'Region': regionName,
                       'Area': '{}'.format(ds.area[index].values),
                       mainRunName: '{}'.format(ds.meltRates[index].values)}
                if dsControl is not None:
                    row[controlRunName] = \
                        '{}'.format(dsControl.meltRates[index].values)
                writer.writerow(row)

        tableFileName = '{}/antarcticMeltFluxTable_{}.csv'.format(outDirectory,
                                                                  self.season)

        with open(tableFileName, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldNames)

            writer.writeheader()
            for index, regionName in enumerate(regionNames):
                row = {'Region': regionName,
                       'Area': '{}'.format(ds.area[index].values),
                       mainRunName: '{}'.format(ds.totalMeltFlux[index].values)}
                if dsControl is not None:
                    row[controlRunName] = \
                        '{}'.format(dsControl.totalMeltFlux[index].values)
                writer.writerow(row)
