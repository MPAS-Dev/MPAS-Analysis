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
import xarray as xr
import datetime

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask


class ClimatologyMapFluxes(AnalysisTask):
    """
    An analysis task for plotting sea surface fluxes. Observational comparison
    is not supported because no observational datasets are currently available.
    """
    # Authors
    # -------
    # Carolyn Begeman

    def __init__(self, config, mpasClimatologyTask, controlConfig=None,
                 fluxType='mass'):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : tranche.Tranche
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : tranche.Tranche, optional
            Configuration options for a control run (if any)

        fluxType : str, optional
            The type of surface fluxes, which corresponds to a different config
            section. One of 'mass' or 'heat'.
        """
        # Authors
        # -------
        # Carolyn Begeman

        taskName = f'climatologyMap{fluxType.title()}Fluxes'
        groupLink = taskName

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapFluxes, self).__init__(
            config=config, taskName=taskName,
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'fluxes', fluxType])

        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getexpression(taskName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {taskName} does not contain '
                             'valid list of seasons')

        comparisonGridNames = config.getexpression(taskName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {taskName} does not contain '
                             'valid list of comparison grids')

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        variableList = config.getexpression(taskName, 'variables')
        for variable in variableList:

            mpasFieldName = f'timeMonthly_avg_{variable}'

            # TemperatureFlux variables have different units and need to be
            # converted to the units of heat flux variables to be compared
            if 'TemperatureFlux' in variable:
                remapClimatologySubtask = RemapMpasTemperatureFluxClimatology(
                    mpasClimatologyTask=mpasClimatologyTask,
                    parentTask=self,
                    climatologyName=variable,
                    variableList=[mpasFieldName],
                    comparisonGridNames=comparisonGridNames,
                    seasons=seasons,
                    iselValues=iselValues,
                    subtaskName=f'remap_{variable}')
                mpasFieldName = 'timeMonthly_avg_' \
                    f'{variable.replace("TemperatureFlux", "HeatFlux")}'
            else:
                remapClimatologySubtask = RemapMpasClimatologySubtask(
                    mpasClimatologyTask=mpasClimatologyTask,
                    parentTask=self,
                    climatologyName=variable,
                    variableList=[mpasFieldName],
                    comparisonGridNames=comparisonGridNames,
                    seasons=seasons,
                    iselValues=iselValues,
                    subtaskName=f'remap_{variable}')

            remapObservationsSubtask = None
            galleryName = variable
            if controlConfig is None:
                refTitleLabel = None
                refFieldName = None
                diffTitleLabel = None
            else:
                control_run_name = controlConfig.get('runs', 'mainRunName')
                refTitleLabel = f'Control: {control_run_name}'
                refFieldName = mpasFieldName
                diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtaskName = f'plot_{variable}_{season}_{comparisonGridName}'
                    outFileName = f'{variable}_{season}_{comparisonGridName}'
                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparisonGridName, remapClimatologySubtask,
                        remapObservationsSubtask, controlConfig=controlConfig,
                        subtaskName=subtaskName)

                    if 'TemperatureFlux' in variable:
                        fieldNameInTitle = variable.replace('TemperatureFlux', 'HeatFlux')
                    else:
                        fieldNameInTitle = variable

                    if 'HeatFlux' in mpasFieldName:
                        groupSubtitle = 'Heat fluxes'
                        unitsLabel = r'W m$^{-2}$'
                    else:
                        groupSubtitle = 'Mass fluxes'
                        unitsLabel = r'kg m$^{-2}$ s^${-1}$'

                    subtask.set_plot_info(
                        outFileLabel=outFileName,
                        fieldNameInTitle=fieldNameInTitle,
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=unitsLabel,
                        imageCaption=variable,
                        galleryGroup='surface fluxes',
                        groupSubtitle=None,
                        groupLink=groupLink,
                        galleryName=galleryName)

                    self.add_subtask(subtask)

# adds to the functionality of RemapMpasClimatology
class RemapMpasTemperatureFluxClimatology(RemapMpasClimatologySubtask):
    """
    A subtask for computing climatologies of heat flux from temperature flux
    """
    # Authors
    # -------
    # Carolyn Begeman

    def customize_masked_climatology(self, climatology, season):
        """
        Construct velocity magnitude as part of the climatology

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
        # Carolyn Begeman

        # first, call the base class's version of this function so we extract
        # the desired slices.
        variable = self.variableList[0]
        climatology = super(RemapMpasTemperatureFluxClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        # calculate heat flux from temperature flux
        scaleFactor = constants.rho_sw * constants.cp_sw  # C m s^-1 to W m^-2
        heatFlux = 0.5 * scaleFactor * climatology[variable]
        # drop unnecessary fields before re-mapping
        climatology.drop_vars([variable])

        # this creates a variable with heat flux units in climatology (like netcdf)
        variable = variable.replace('TemperatureFlux', 'HeatFlux')
        climatology[variable] = heatFlux
        climatology[variable].attrs['units'] = 'W m$^[-2]$'
        climatology[variable].attrs['description'] = \
            f'{variable} converted to heat flux'

        return climatology
