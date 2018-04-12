from __future__ import absolute_import, division, print_function, \
    unicode_literals

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.ocean.compute_transects_subtask import \
    ComputeTransectsSubtask, TransectsObservations

from mpas_analysis.ocean.plot_transect_subtask import PlotTransectSubtask

from mpas_analysis.shared.io.utility import build_config_full_path

from collections import OrderedDict


class WoceTransects(AnalysisTask):  # {{{
    """
    Plot model output at WOCE transects and compare it against WOCE
    observations

    Attributes
    ----------

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, mpasClimatologyTask, refConfig=None):

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

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)

        Authors
        -------
        Xylar Asay-Davis
        '''

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

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'woceSubdirectory')

        obsFileNames = OrderedDict()
        for transectName in ['WOCE_A21_Drake_Passage',
                             'WOCE_A23_South_Atlantic',
                             'WOCE_A12_Prime_Meridian']:

            fileName = '{}/{}.nc'.format(observationsDirectory, transectName)
            obsFileNames[transectName] = fileName

        fields = \
            {'temperature':
                {'mpas': 'timeMonthly_avg_activeTracers_temperature',
                 'obs': 'potentialTemperature',
                 'units': r'$^\circ$C'},
             'salinity':
                {'mpas': 'timeMonthly_avg_activeTracers_salinity',
                 'obs': 'salinity',
                 'units': r'PSU'}}

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

        plotObs = refConfig is None
        if plotObs:

            refTitleLabel = 'Observations (WOCE)'

            galleryName = 'Observations: AVISO'

            diffTitleLabel = 'Model - Observations'

        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Ref: {}'.format(refRunName)

            diffTitleLabel = 'Main - Reference'

        for transectName in obsFileNames:
            for season in seasons:
                for fieldName in fields:
                    outFileLabel = fieldName
                    if plotObs:
                        refFieldName = fields[fieldName]['obs']
                    else:
                        refFieldName = fields[fieldName]['mpas']

                    fieldNameUpper = fieldName[0].upper() + fieldName[1:]
                    fieldNameInTytle = '{} from {}'.format(
                            fieldNameUpper,
                            transectName.replace('_',' '))
                    galleryGroup = 'WOCE {}'.format(fieldNameUpper)
                    groupLink = 'woce_{}'.format(fieldName)

                    # make a new subtask for this season and comparison grid
                    subtask = PlotTransectSubtask(self, season, transectName,
                                                  fieldName,
                                                  computeTransectsSubtask,
                                                  plotObs, refConfig)

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
                            galleryGroup=galleryGroup,
                            groupSubtitle=None,
                            groupLink=groupLink,
                            galleryName=galleryName,
                            configSectionName='woce{}Transects'.format(
                                    fieldNameUpper))

                    self.add_subtask(subtask)
        # }}}

    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
