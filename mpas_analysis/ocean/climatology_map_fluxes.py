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

from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask


class ClimatologyMapFluxes(AnalysisTask):
    """
    An analysis task for plotting sea surface fluxes
    """
    # Authors
    # -------
    # Carolyn Begeman

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
        # Authors
        # -------
        # Carolyn Begeman

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapFluxes, self).__init__(
            config=config, taskName=f'climatologyMapFluxes',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', 'fluxes'])

        sectionName = self.taskName

        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError(f'config section {sectionName} does not contain valid list '
                             'of seasons')

        comparisonGridNames = config.getexpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError(f'config section {sectionName} does not contain valid list '
                             'of comparison grids')

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        variableList = config.getexpression(sectionName, 'variables')
        unitsLabels = {}
        mpasVariableList = []
        remapClimatologySubtasks = {}
        groupLink = 'ocnfluxclimo'
        for variable in variableList:

            unitsLabels[variable] = r'kg m$^2$ s^${-1}$'
            mpasFieldName = f'timeMonthly_avg_{variable}'
            mpasVariableList.append(mpasFieldName)

            remapClimatologySubtasks[variable] = RemapMpasClimatologySubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName=variable,
                variableList=[mpasFieldName],
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                iselValues=iselValues,
                subtaskName=f'remap_{variable}')

            remapObservationsSubtask = None
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
                galleryName = variable
                unitsLabel = unitsLabels[variable]
                for season in seasons:
                    # make a new subtask for this season and comparison grid
                    subtaskName = f'plot_{variable}_{season}_{comparisonGridName}'
                    outFileName = f'{variable}_{season}_{comparisonGridName}'
                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparisonGridName, remapClimatologySubtasks[variable],
                        remapObservationsSubtask, controlConfig=controlConfig,
                        subtaskName=subtaskName)

                    subtask.set_plot_info(
                        outFileLabel=outFileName,
                        fieldNameInTitle=variable,
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=unitsLabel,
                        imageCaption=variable,
                        galleryGroup='surface fluxes',
                        groupSubtitle='Mass fluxes',
                        groupLink=groupLink,
                        galleryName=galleryName)

                    self.add_subtask(subtask)
