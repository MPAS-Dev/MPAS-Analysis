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
from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_obs_path

from collections import OrderedDict


class WoceTransects(AnalysisTask):  # {{{
    """
    Plot model output at WOCE transects and compare it against WOCE
    observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):

        # {{{
        '''
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        controlConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a control run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        tags = ['climatology', 'transect', 'woce']

        # call the constructor from the base class (AnalysisTask)
        super(WoceTransects, self).__init__(
            config=config, taskName='woceTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        seasons = config.getExpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getExpression(
                sectionName, 'verticalComparisonGrid', usenumpyfunc=True)

        observationsDirectory = build_obs_path(
            config, 'ocean', 'woceSubdirectory')

        obsFileNames = OrderedDict()

        obsFileNames['WOCE_A21_Drake_Passage'] = \
            'WOCE_A21_Drake_Passage_20181126.nc'
        obsFileNames['WOCE_A23_South_Atlantic'] = \
            'WOCE_A23_South_Atlantic_20181126.nc'
        obsFileNames['WOCE_A12_Prime_Meridian'] = \
            'WOCE_A12_Prime_Meridian_20181126.nc'

        for transectName in obsFileNames:
            fileName = '{}/{}'.format(observationsDirectory,
                                      obsFileNames[transectName])
            obsFileNames[transectName] = fileName

        fields = \
            {'temperature':
                {'mpas': 'timeMonthly_avg_activeTracers_temperature',
                 'obs': 'potentialTemperature',
                 'titleName': 'Potential Temperature',
                 'units': r'$\degree$C'},
             'salinity':
                {'mpas': 'timeMonthly_avg_activeTracers_salinity',
                 'obs': 'salinity',
                 'titleName': 'Salinity',
                 'units': r'PSU'},
             'potentialDensity':
                {'mpas': 'timeMonthly_avg_potentialDensity',
                 'obs': 'potentialDensity',
                 'titleName': 'Potential Density',
                 'units': r'kg m$^{-3}$'}}

        transectCollectionName = 'WOCE_transects'
        if horizontalResolution != 'obs':
            transectCollectionName = '{}_{}km'.format(transectCollectionName,
                                                      horizontalResolution)

        transectsObservations = TransectsObservations(config, obsFileNames,
                                                      horizontalResolution,
                                                      transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='WOCE',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = controlConfig is None
        if plotObs:

            refTitleLabel = 'Observations (WOCE)'

            diffTitleLabel = 'Model - Observations'

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = 'Control: {}'.format(controlRunName)

            diffTitleLabel = 'Main - Control'

        fieldNameDict = {'temperature': 'temperatureTransect',
                         'salinity': 'salinityTransect',
                         'potentialDensity': 'potentialDensityTransect'}

        for fieldName in fields:
            for transectName in obsFileNames:
                for season in seasons:
                    outFileLabel = fieldNameDict[fieldName]
                    if plotObs:
                        refFieldName = fields[fieldName]['obs']
                    else:
                        refFieldName = fields[fieldName]['mpas']

                    fieldNameUpper = fieldName[0].upper() + fieldName[1:]
                    titleName = fields[fieldName]['titleName']
                    fieldNameInTytle = '{} from {}'.format(
                        titleName,
                        transectName.replace('_', ' '))

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(self, season, transectName,
                                                  fieldName,
                                                  computeTransectsSubtask,
                                                  plotObs, controlConfig)

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTytle,
                        mpasFieldName=fields[fieldName]['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=fields[fieldName]['units'],
                        imageCaption='{} {}'.format(fieldNameInTytle,
                                                    season),
                        galleryGroup='WOCE Transects',
                        groupSubtitle=None,
                        groupLink='woce',
                        galleryName=titleName,
                        configSectionName='woce{}Transects'.format(
                            fieldNameUpper))

                    self.add_subtask(subtask)
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
