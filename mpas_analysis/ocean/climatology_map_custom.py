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

import numpy as np
import xarray as xr
from mpas_tools.cime.constants import constants as cime_constants

from mpas_analysis.ocean.remap_depth_slices_subtask import (
    RemapDepthSlicesSubtask
)
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask


class ClimatologyMapCustom(AnalysisTask):
    """
    A flexible analysis task for plotting climatologies of any MPAS-Ocean field
    on cells from timeSeriesStatsMonthly at various depths (if the field has
    vertical levels) and for various seasons.

    Various derived fields are also supported:

    * velocity magnitude
    * thermal forcing (temperature - freezing temperature)
    """

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """

        taskName = 'climatologyMapCustom'

        sectionName = taskName
        variablesNames = config.getexpression(sectionName, 'variables')

        tags = ['climatology', 'horizontalMap'] + variablesNames

        # call the constructor from the base class (AnalysisTask)
        super().__init__(
            config=config, taskName=taskName, componentName='ocean', tags=tags)

        if len(variablesNames) == 0:
            return

        variablesNames = config.getexpression(sectionName, 'variables')
        if len(variablesNames) == 0:
            return

        availableVariables = config.getexpression(sectionName,
                                                  'availableVariables')

        variables = {varName: availableVariables[varName] for varName in
                     variablesNames}

        for varName in variablesNames:
            if 'has_depth' not in variables[varName]:
                # we assume variables have depth unless otherwise specified
                variables[varName]['has_depth'] = True

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             'valid list of seasons')

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {sectionName} does not contain '
                             'valid list of comparison grids')

        depths = config.getexpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError(f'config section {sectionName} does not '
                             f'contain valid list of depths')

        remapMpasSubtask = RemapMpasDerivedVariableClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='custom3D',
            variables=variables,
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames)

        galleryGroup = 'Custom Climatology Maps'
        groupLink = 'custclimmaps'

        for varName, metadata in variables.items():
            title = metadata['title']
            units = metadata['units']
            hasDepth = metadata['has_depth']
            upperVarName = varName[0].upper() + varName[1:]
            varSectionName = f'{self.taskName}{upperVarName}'

            remapObsSubtask = None

            refTitleLabel = None
            diffTitleLabel = None
            if controlConfig is not None:
                controlRunName = controlConfig.get('runs', 'mainRunName')
                refTitleLabel = f'Control: {controlRunName}'
                diffTitleLabel = 'Main - Control'

            if hasDepth:
                localDepths = depths
            else:
                localDepths = [None]

            for comparisonGridName in comparisonGridNames:
                for depth in localDepths:
                    for season in seasons:

                        subtaskName = f'plot{upperVarName}_{season}_' \
                                      f'{comparisonGridName}'
                        if depth is not None:
                            subtaskName = f'{subtaskName}_depth_{depth}'

                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapMpasSubtask,
                            remapObsClimatologySubtask=remapObsSubtask,
                            controlConfig=controlConfig,
                            depth=depth,
                            subtaskName=subtaskName)

                        subtask.set_plot_info(
                            outFileLabel=f'cust_{varName}',
                            fieldNameInTitle=title,
                            mpasFieldName=varName,
                            refFieldName=varName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=units,
                            imageCaption=title,
                            galleryGroup=galleryGroup,
                            groupSubtitle=None,
                            groupLink=groupLink,
                            galleryName=title,
                            configSectionName=varSectionName)

                        self.add_subtask(subtask)


class RemapMpasDerivedVariableClimatology(RemapDepthSlicesSubtask):
    """
    A subtask for computing derived variables (such as velocity magnitude and
    thermal forcing) as part of remapping climatologies at depth slices

    Attributes
    ----------
    variables : dict of dict
        A dictionary of variable definitions, with variable names as keys

    """

    def __init__(self, mpasClimatologyTask, parentTask, climatologyName,
                 variables, seasons, depths, comparisonGridNames):

        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        mpasClimatologyTask : MpasClimatologyTask
            The task that produced the climatology to be remapped

        parentTask : AnalysisTask
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        climatologyName : str
            A name that describes the climatology (e.g. a short version of
            the important field(s) in the climatology) used to name the
            subdirectories for each stage of the climatology

        variables : dict of dict
            A dictionary of variable definitions, with variable names as keys

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        depths : list of {None, float, 'top', 'bot'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.
        """
        self.variables = variables

        mpasVariables = set()
        for variable in variables.values():
            for mpasVariable in variable['mpas']:
                mpasVariables.add(mpasVariable)

        # call the constructor from the base class
        # (RemapMpasClimatologySubtask)
        super().__init__(
            mpasClimatologyTask, parentTask, climatologyName, mpasVariables,
            seasons, depths, comparisonGridNames, iselValues=None)

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

        # first, compute the derived variables, which may rely on having the
        # full 3D variables available

        derivedVars = []
        self._add_vel_mag(climatology, derivedVars)
        self._add_thermal_forcing(climatology, derivedVars)

        # then, call the superclass's version of this function so we extract
        # the desired slices (but before renaming because it expects the
        # original MPAS variable names)
        climatology = super().customize_masked_climatology(climatology,
                                                           season)
        # finally, rename the variables and add metadata
        for varName, variable in self.variables.items():
            if varName not in derivedVars:
                # rename variables from MPAS names to shorter names
                mpasvarNames = variable['mpas']
                if len(mpasvarNames) == 1:
                    mpasvarName = mpasvarNames[0]
                    climatology[varName] = climatology[mpasvarName]
                    climatology.drop_vars(mpasvarName)

            climatology[varName].attrs['units'] = variable['units']
            climatology[varName].attrs['description'] = variable['title']

        return climatology

    def _add_vel_mag(self, climatology, derivedVars):
        """
        Add the velocity magnitude to the climatology if requested
        """
        varName = 'velocityMagnitude'
        if varName not in self.variables:
            return

        derivedVars.append(varName)

        zonalVel = climatology.timeMonthly_avg_velocityZonal
        meridVel = climatology.timeMonthly_avg_velocityMeridional
        climatology[varName] = np.sqrt(zonalVel**2 + meridVel**2)

    def _add_thermal_forcing(self, climatology, derivedVars):
        """
        Add thermal forcing to the climatology if requested
        """
        varName = 'thermalForcing'
        if varName not in self.variables:
            return

        derivedVars.append(varName)

        c0 = self.namelist.getfloat(
            'config_land_ice_cavity_freezing_temperature_coeff_0')
        cs = self.namelist.getfloat(
            'config_land_ice_cavity_freezing_temperature_coeff_S')
        cp = self.namelist.getfloat(
            'config_land_ice_cavity_freezing_temperature_coeff_p')
        cps = self.namelist.getfloat(
            'config_land_ice_cavity_freezing_temperature_coeff_pS')

        temp = climatology.timeMonthly_avg_activeTracers_temperature
        salin = climatology.timeMonthly_avg_activeTracers_salinity
        dens = climatology.timeMonthly_avg_density
        thick = climatology.timeMonthly_avg_layerThickness

        dp = cime_constants['SHR_CONST_G']*dens*thick
        press = dp.cumsum(dim='nVertLevels') - 0.5*dp

        # add land ice pressure if available
        ds_restart = xr.open_dataset(self.restartFileName)
        ds_restart = ds_restart.isel(Time=0)
        if 'landIcePressure' in ds_restart:
            press += ds_restart.landIcePressure

        tempFreeze = c0 + cs*salin + cp*press + cps*press*salin

        climatology[varName] = temp - tempFreeze
