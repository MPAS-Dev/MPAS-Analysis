# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class OsnapTransects(AnalysisTask):
    """
    Plot model output at OSNAP transects and compare it against OSNAP
    observations
    """
    # Authors
    # -------
    # Alice Barthel, Xylar Asay-Davis

    def __init__(self, config, mpasClimatologyTask, controlConfig=None):
        """
        Construct the analysis task and adds it as a subtask of the
        ``parentTask``.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted
            as a transect

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Alice Barthel

        tags = ['climatology', 'transect', 'osnap', 'publicObs', 'arctic']

        # call the constructor from the base class (AnalysisTask)
        super(OsnapTransects, self).__init__(
            config=config, taskName='osnapTransects',
            componentName='ocean',
            tags=tags)

        sectionName = self.taskName

        seasons = config.getexpression(sectionName, 'seasons')

        horizontalResolution = config.get(sectionName, 'horizontalResolution')

        verticalComparisonGridName = config.get(sectionName,
                                                'verticalComparisonGridName')

        if verticalComparisonGridName in ['mpas', 'obs']:
            verticalComparisonGrid = None
        else:
            verticalComparisonGrid = config.getexpression(
                sectionName, 'verticalComparisonGrid', use_numpyfunc=True)

        verticalBounds = config.getexpression(sectionName, 'verticalBounds')

        horizontalBounds = config.getexpression(
            sectionName, 'horizontalBounds')

        observationsDirectory = build_obs_path(
            config, 'ocean', 'osnapSubdirectory')

        origObsFileNames = \
            {'OSNAP_East': 'OSNAP_East_ANN_20220926.nc',
             'OSNAP_West': 'OSNAP_West_ANN_20220926.nc'}

        obsFileNames = {}
        for transectName in horizontalBounds:
            found = False
            for name in origObsFileNames:
                if transectName.startswith(name):
                    obsFileNames[transectName] = origObsFileNames[name]
                    found = True
                    break
            if not found:
                raise ValueError(f'Keys for horizontalBounds must start '
                                 f'with one of {list(origObsFileNames)}')

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
                 'units': r'PSU'}}

        transectCollectionName = 'OSNAP_transects'
        if horizontalResolution not in ['obs', 'mpas']:
            transectCollectionName = \
                f'{transectCollectionName}_{horizontalResolution}km'

        transectsObservations = TransectsObservations(config, obsFileNames,
                                                      horizontalResolution,
                                                      transectCollectionName)

        computeTransectsSubtask = ComputeTransectsSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='OSNAP',
            transectCollectionName=transectCollectionName,
            variableList=[field['mpas'] for field in fields.values()],
            seasons=seasons,
            obsDatasets=transectsObservations,
            verticalComparisonGridName=verticalComparisonGridName,
            verticalComparisonGrid=verticalComparisonGrid)

        plotObs = controlConfig is None
        if plotObs:

            refTitleLabel = 'Observations (OSNAP)'

            diffTitleLabel = 'Model - Observations'

        else:
            controlRunName = controlConfig.get('runs', 'mainRunName')
            refTitleLabel = f'Control: {controlRunName}'

            diffTitleLabel = 'Main - Control'

        fieldNameDict = {'temperature': 'temperatureTransect',
                         'salinity': 'salinityTransect'}

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
                    fieldNameInTitle = \
                        f'{titleName} from {transectName.replace("_", " ")}'

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(
                        self, season, transectName, fieldName,
                        computeTransectsSubtask, plotObs, controlConfig,
                        horizontalBounds[transectName])

                    subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle=fieldNameInTitle,
                        mpasFieldName=fields[fieldName]['mpas'],
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=fields[fieldName]['units'],
                        imageCaption=f'{fieldNameInTitle} {season}',
                        galleryGroup='OSNAP Transects',
                        groupSubtitle=None,
                        groupLink='osnap',
                        galleryName=titleName,
                        configSectionName=f'osnap{fieldNameUpper}Transects',
                        verticalBounds=verticalBounds)

                    self.add_subtask(subtask)
