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
from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask


class ClimatologyMapMLDMinMax(AnalysisTask):
    """
    An analysis task for comparison of mixed layer depth (mld) against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasClimatologyTasks, controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTasks : dict of ``MpasClimatologyTask``
            The tasks that produced the climatology of monthly min and max to
            be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'mld'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapMLDMinMax, self).__init__(
            config=config, taskName='climatologyMapMLDMinMax',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName, 'publicObs',
                  'min', 'max'])

        self._add_tasks(config, mpasClimatologyTasks, controlConfig,
                        title='Density MLD',
                        mpasVariableSuffix='dThreshMLD',
                        filePrefix='d_mld',
                        sectionPrefix='climatologyMapDensityMLD')
        self._add_tasks(config, mpasClimatologyTasks, controlConfig,
                        title='Temperature MLD',
                        mpasVariableSuffix='tThreshMLD',
                        filePrefix='t_mld',
                        sectionPrefix='climatologyMapTemperatureMLD')
        self._add_tasks(config, mpasClimatologyTasks, controlConfig,
                        title='Boundary-Layer Depth',
                        mpasVariableSuffix='boundaryLayerDepth',
                        filePrefix='bld',
                        sectionPrefix='climatologyMapBLD')

    def setup_and_check(self):
        """
        Check if MLD capability was turned on in the run.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapMLDMinMax, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_mixedlayerdepths_enable',
            raiseException=True)

        self.check_analysis_enabled(
            analysisOptionName='config_AM_timeSeriesStatsMonthlyMin_enable',
            raiseException=True)
        self.check_analysis_enabled(
            analysisOptionName='config_AM_timeSeriesStatsMonthlyMax_enable',
            raiseException=True)

    def _add_tasks(self, config, mpasClimatologyTasks, controlConfig,
                   title, mpasVariableSuffix, filePrefix, sectionPrefix):
        """
        Add tasks for a given variable
        """
        iselValues = None

        # read in what seasons we want to plot
        seasons = config.getexpression(self.taskName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(self.taskName))

        comparisonGridNames = config.getexpression(self.taskName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(self.taskName))

        remapClimatologySubtasks = {}
        mpasFieldNames = {}
        for op in ['min', 'max']:
            upperOp = op[0].upper() + op[1:]

            mpasFieldNames[op] = 'timeMonthly{}_{}_{}'.format(
                upperOp, op, mpasVariableSuffix)
            # the variable 'timeMonthly_avg_dThreshMLD' will be added to
            # mpasClimatologyTask along with the seasons.
            remapClimatologySubtasks[op] = RemapMpasClimatologySubtask(
                mpasClimatologyTask=mpasClimatologyTasks[op],
                parentTask=self,
                climatologyName=filePrefix,
                variableList=[mpasFieldNames[op]],
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                iselValues=iselValues,
                subtaskName='remap{}_{}'.format(upperOp, filePrefix))

        if controlConfig is None:
            fieldNameInTitle = 'Max {}'.format(title)
            refTitleLabel = 'Min {}'.format(title)
            diffTitleLabel = 'Max - Min {}'.format(title)
            galleryName = title
            outFileLabel = '{}_min_max'.format(filePrefix)
            sectionName = '{}MinMax'.format(sectionPrefix)

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    subtaskName = 'plot_{}_{}_{}'.format(filePrefix,
                                                         season,
                                                         comparisonGridName)
                    # make a new subtask for this season and comparison grid
                    firstTask = remapClimatologySubtasks['max']
                    secondTask = remapClimatologySubtasks['min']
                    subtask = PlotClimatologyMapSubtask(
                        self, season, comparisonGridName,
                        remapMpasClimatologySubtask=firstTask,
                        secondRemapMpasClimatologySubtask=secondTask,
                        controlConfig=controlConfig,
                        subtaskName=subtaskName)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTitle,
                        mpasFieldName=mpasFieldNames['max'],
                        refFieldName=mpasFieldNames['min'],
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'm',
                        imageCaption='Min/Max Mixed-Layer Depth',
                        galleryGroup='Min/Max Mixed-Layer Depth',
                        groupSubtitle=None,
                        groupLink='mld_min_max',
                        galleryName=galleryName,
                        configSectionName=sectionName)

                    self.add_subtask(subtask)

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = title
            refTitleLabel = 'Control: {}'.format(controlRunName)

            diffTitleLabel = 'Main - Control'

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    for op in ['min', 'max']:
                        upperOp = op[0].upper() + op[1:]
                        subtaskName = 'plot{}_{}_{}_{}'.format(
                            upperOp, filePrefix, season, comparisonGridName)
                        fieldNameInTitle = '{} {}'.format(upperOp, title)
                        outFileLabel = '{}_{}'.format(filePrefix, op)
                        sectionName = '{}{}'.format(sectionPrefix, upperOp)

                        # make a new subtask for this season and comparison
                        # grid
                        subtask = PlotClimatologyMapSubtask(
                            self, season, comparisonGridName,
                            remapClimatologySubtasks[op],
                            controlConfig=controlConfig,
                            subtaskName=subtaskName)

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=fieldNameInTitle,
                            mpasFieldName=mpasFieldNames[op],
                            refFieldName=mpasFieldNames[op],
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=r'm',
                            imageCaption=fieldNameInTitle,
                            galleryGroup='Min/Max Mixed-Layer Depth',
                            groupSubtitle=None,
                            groupLink='mld_{}'.format(op),
                            galleryName=galleryName,
                            configSectionName=sectionName)

                        self.add_subtask(subtask)
