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
"""
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.
"""
# Authors
# -------
# Carolyn Begeman

import numpy

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.shared.plot import PlotClimatologyMapSubtask
from mpas_analysis.ocean.remap_sose_climatology import RemapSoseClimatology

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapVel(AnalysisTask):
    """
    An analysis task for comparison of antarctic field against the Southern
    Ocean State Estimate
    """
    # Authors
    # -------
    # Carolyn Begeman

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):
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

        fields = \
            [{'prefix': 'zonalVelocity',
              'mpas': 'timeMonthly_avg_velocityZonal',
              'units': r'm s$^{-1}$',
              'titleName': 'Zonal Velocity',
              '3D': True,
              'obsFilePrefix': 'zonal_vel',
              'obsFieldName': 'zonalVel',
              'obsBotFieldName': 'botZonalVel'},
             {'prefix': 'meridionalVelocity',
              'mpas': 'timeMonthly_avg_velocityMeridional',
              'units': r'm s$^{-1}$',
              'titleName': 'Meridional Velocity',
              '3D': True,
              'obsFilePrefix': 'merid_vel',
              'obsFieldName': 'meridVel',
              'obsBotFieldName': 'botMeridVel'},
             {'prefix': 'velocityMagnitude',
              'mpas': 'velMag',
              'units': r'm s$^{-1}$',
              'titleName': 'Velocity Magnitude',
              '3D': True,
              'obsFilePrefix': 'vel_mag',
              'obsFieldName': 'velMag',
              'obsBotFieldName': 'botVelMag'}]

        tags = ['climatology', 'horizontalMap', 'climatologyMapVel']
        for field in fields:
            tags.append(field['prefix'])

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapVel, self).__init__(
            config=config, taskName='climatologyMapVel',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        fileSuffix = config.get(sectionName, 'fileSuffix')
        if fileSuffix.endswith('.nc'):
            fileSuffix = fileSuffix.strip('.nc')

        fieldList = config.getexpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

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

        if not numpy.any([field['3D'] for field in fields]):
            depths = None
        else:
            depths = config.getexpression(sectionName, 'depths')

            if len(depths) == 0:
                raise ValueError(f'config section {sectionName} does not '
                                 'contain valid list of depths')

        variableList = []
        for field in fields:
            if field['mpas'] != 'velMag':
                variableList.append(field['mpas'])
            else:
                for variable_dependency in ['timeMonthly_avg_velocityZonal',
                                            'timeMonthly_avg_velocityMeridional']:
                    if variable_dependency not in variableList:
                        variableList.append(variable_dependency)

        shallowVsDeepColormapDepth = config.getfloat(
            sectionName, 'shallowVsDeepColormapDepth')

        shallow = []
        for depth in depths:
            if depth == 'top':
                shallow.append(True)
            elif depth == 'bot':
                shallow.append(False)
            else:
                shallow.append(depth >= shallowVsDeepColormapDepth)

        if depths is None:
            remapMpasSubtask = RemapMpasClimatologySubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='vel',
                variableList=variableList,
                seasons=seasons,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)
        else:
            remapMpasSubtask = RemapMpasVelMagClimatology(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='vel',
                variableList=variableList,
                seasons=seasons,
                depths=depths,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            sectionName = f'{self.taskName}{upperFieldPrefix}'

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = None

            remapObsSubtask = None

            refFieldName = None
            outFileLabel = fieldPrefix
            refTitleLabel = None
            diffTitleLabel = None
            galleryName = 'Velocity'
            if controlConfig is not None:
                refFieldName = field['mpas']
                controlRunName = controlConfig.get('runs', 'mainRunName')
                galleryName = f'Control: {controlRunName}'
                refTitleLabel = galleryName
                diffTitleLabel = 'Main - Control'

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = [None]

            for comparisonGridName in comparisonGridNames:
                for depthIndex, depth in enumerate(fieldDepths):
                    for season in seasons:

                        subtaskName = f'plot{upperFieldPrefix}_{season}_' \
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

                        configSectionName = 'climatologyMapVel_' \
                                            f'{upperFieldPrefix}'

                        # if available, use a separate color map for shallow
                        # and deep
                        if depth is not None:
                            if shallow[depthIndex]:
                                suffix = 'Shallow'
                            else:
                                suffix = 'Deep'
                            testSectionName = f'{configSectionName}{suffix}'

                        if config.has_section(testSectionName):
                            configSectionName = testSectionName
                        else:
                            configSectionName = 'climatologyMapVel'

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=field['titleName'],
                            mpasFieldName=field['mpas'],
                            refFieldName=refFieldName,
                            refTitleLabel=refTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=field['units'],
                            imageCaption=field['titleName'],
                            galleryGroup=galleryName,
                            groupSubtitle=None,
                            groupLink=f'{fieldPrefix}Vel',
                            galleryName=field['titleName'],
                            configSectionName=configSectionName)

                        self.add_subtask(subtask)


class RemapMpasVelMagClimatology(RemapDepthSlicesSubtask):
    """
    A subtask for computing climatologies of velocity magnitude from zonal
    and meridional components
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
        climatology = super(RemapMpasVelMagClimatology,
                            self).customize_masked_climatology(climatology,
                                                               season)

        if 'timeMonthly_avg_velocityZonal' in climatology and \
                'timeMonthly_avg_velocityMeridional' in climatology:
            zonalVel = climatology.timeMonthly_avg_velocityZonal
            meridVel = climatology.timeMonthly_avg_velocityMeridional
            climatology['velMag'] = numpy.sqrt(zonalVel**2 + meridVel**2)
            climatology.velMag.attrs['units'] = 'm s$^{-1}$'
            climatology.velMag.attrs['description'] = 'velocity magnitude'

        return climatology
