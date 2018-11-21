# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.
'''
# Authors
# -------
# Xylar Asay-Davis

import numpy

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask

from mpas_analysis.ocean.remap_depth_slices_subtask import \
    RemapDepthSlicesSubtask
from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask
from mpas_analysis.ocean.remap_sose_climatology import RemapSoseClimatology

from mpas_analysis.shared.io.utility import build_config_full_path


class ClimatologyMapSose(AnalysisTask):  # {{{
    """
    An analysis task for comparison of antarctic field against the Southern
    Ocean State Estimate
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask,
                 refConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fields = \
            [{'prefix': 'temperature',
              'mpas': 'timeMonthly_avg_activeTracers_temperature',
              'units': r'$\degree$C',
              'titleName': 'Potential Temperature',
              '3D': True,
              'obsFilePrefix': 'pot_temp',
              'obsFieldName': 'theta',
              'obsBotFieldName': 'botTheta'},
             {'prefix': 'salinity',
              'mpas': 'timeMonthly_avg_activeTracers_salinity',
              'units': r'PSU',
              'titleName': 'Salinity',
              '3D': True,
              'obsFilePrefix': 'salinity',
              'obsFieldName': 'salinity',
              'obsBotFieldName': 'botSalinity'},
             {'prefix': 'potentialDensity',
              'mpas': 'timeMonthly_avg_potentialDensity',
              'units': r'kg m$^{-3}$',
              'titleName': 'Potential Density',
              '3D': True,
              'obsFilePrefix': 'pot_den',
              'obsFieldName': 'potentialDensity',
              'obsBotFieldName': 'botPotentialDensity'},
             {'prefix': 'mixedLayerDepth',
              'mpas': 'timeMonthly_avg_dThreshMLD',
              'units': r'm',
              'titleName': 'Mixed Layer Depth',
              '3D': False,
              'obsFilePrefix': 'mld',
              'obsFieldName': 'mld',
              'obsBotFieldName': None},
             {'prefix': 'zonalVelocity',
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

        tags = ['climatology', 'horizontalMap', 'sose', 'publicObs'] \
            + [field['prefix'] for field in fields]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSose, self).__init__(
                config=config, taskName='climatologyMapSose',
                componentName='ocean',
                tags=tags)

        sectionName = self.taskName

        fieldList = config.getExpression(sectionName, 'fieldList')
        fields = [field for field in fields if field['prefix'] in fieldList]

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of comparison grids'.format(
                                     sectionName))

        if not numpy.any([field['3D'] for field in fields]):
            depths = None
        else:
            depths = config.getExpression(sectionName, 'depths')

            if len(depths) == 0:
                raise ValueError('config section {} does not contain valid '
                                 'list of depths'.format(sectionName))

        variableList = [field['mpas'] for field in fields
                        if field['mpas'] != 'velMag']

        if depths is None:
            remapMpasSubtask = RemapMpasClimatologySubtask(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='SOSE',
                variableList=variableList,
                seasons=seasons,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)
        else:
            remapMpasSubtask = RemapMpasVelMagClimatology(
                mpasClimatologyTask=mpasClimatologyTask,
                parentTask=self,
                climatologyName='SOSE',
                variableList=variableList,
                seasons=seasons,
                depths=depths,
                comparisonGridNames=comparisonGridNames,
                iselValues=None)

        for field in fields:
            fieldPrefix = field['prefix']
            upperFieldPrefix = fieldPrefix[0].upper() + fieldPrefix[1:]
            sectionName = '{}{}'.format(self.taskName, upperFieldPrefix)

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = None

            if refConfig is None:

                controlTitleLabel = 'State Estimate (SOSE)'

                observationsDirectory = build_config_full_path(
                    config, 'oceanObservations', 'soseSubdirectory')

                obsFileName = \
                    '{}/SOSE_2005-2010_monthly_{}_6000.0x' \
                    '6000.0km_10.0km_Antarctic_stereo.nc'.format(
                            observationsDirectory, field['obsFilePrefix'])
                controlFieldName = field['obsFieldName']
                outFileLabel = '{}SOSE'.format(fieldPrefix)
                galleryName = 'State Estimate: SOSE'
                diffTitleLabel = 'Model - State Estimate'

                remapObsSubtask = RemapSoseClimatology(
                        parentTask=self, seasons=seasons, fileName=obsFileName,
                        outFilePrefix='{}SOSE'.format(controlFieldName),
                        fieldName=controlFieldName,
                        botFieldName=field['obsBotFieldName'],
                        depths=fieldDepths,
                        comparisonGridNames=comparisonGridNames,
                        subtaskName='remapObservations{}'.format(
                                upperFieldPrefix))

                self.add_subtask(remapObsSubtask)

            else:
                remapObsSubtask = None
                refRunName = refConfig.get('runs', 'mainRunName')
                galleryName = 'Ref: {}'.format(refRunName)
                controlTitleLabel = galleryName

                controlFieldName = field['mpas']
                outFileLabel = fieldPrefix
                diffTitleLabel = 'Main - Reference'

            if field['3D']:
                fieldDepths = depths
            else:
                fieldDepths = [None]

            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    for depth in fieldDepths:

                        subtaskName = 'plot{}_{}_{}'.format(upperFieldPrefix,
                                                            season,
                                                            comparisonGridName)
                        if depth is not None:
                            subtaskName = '{}_depth_{}'.format(subtaskName,
                                                               depth)

                        subtask = PlotClimatologyMapSubtask(
                            parentTask=self,
                            season=season,
                            comparisonGridName=comparisonGridName,
                            remapMpasClimatologySubtask=remapMpasSubtask,
                            remapObsClimatologySubtask=remapObsSubtask,
                            refConfig=refConfig,
                            depth=depth,
                            subtaskName=subtaskName)

                        subtask.set_plot_info(
                            outFileLabel=outFileLabel,
                            fieldNameInTitle=field['titleName'],
                            mpasFieldName=field['mpas'],
                            controlFieldName=controlFieldName,
                            controlTitleLabel=controlTitleLabel,
                            diffTitleLabel=diffTitleLabel,
                            unitsLabel=field['units'],
                            imageCaption=field['titleName'],
                            galleryGroup=field['titleName'],
                            groupSubtitle=None,
                            groupLink='{}Sose'.format(fieldPrefix),
                            galleryName=galleryName,
                            configSectionName='climatologyMapSose{}'.format(
                                upperFieldPrefix))

                        self.add_subtask(subtask)
        # }}}

    # }}}


class RemapMpasVelMagClimatology(RemapDepthSlicesSubtask):  # {{{
    """
    A subtask for computing climatologies of velocity magnitude from zonal
    and meridional components
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def customize_masked_climatology(self, climatology, season):  # {{{
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
        # Xylar Asay-Davis

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

        return climatology  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
