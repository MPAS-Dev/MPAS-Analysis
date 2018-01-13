
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from ..shared import AnalysisTask

from ..shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from .plot_climatology_map_subtask import PlotClimatologyMapSubtask

from ..shared.io.utility import build_config_full_path

from ..shared.grid import LatLonGridDescriptor


class ClimatologyMapSeaIceThick(AnalysisTask):  # {{{
    """
    An analysis task for comparison of sea ice thickness against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, mpasClimatologyTask, hemisphere,
                 refConfig=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)

        Authors
        -------
        Xylar Asay-Davis
        """
        taskName = 'climatologyMapSeaIceThick{}'.format(hemisphere)

        fieldName = 'seaIceThick'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(
                config=config, taskName=taskName,
                componentName='seaIce',
                tags=['climatology', 'horizontalMap', fieldName])

        mpasFieldName = 'timeMonthly_avg_iceVolumeCell'
        iselValues = None

        sectionName = taskName

        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getExpression(sectionName,
                                                   'comparisonGrids')

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # the variable self.mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName='{}{}'.format(fieldName, hemisphere),
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if refConfig is None:
            refTitleLabel = 'Observations (ICESat)'
            galleryName = 'Observations: ICESat'
            diffTitleLabel = 'Model - Observations'
            refFieldName = 'seaIceThick'
        else:
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Ref: {}'.format(refRunName)
            refFieldName = mpasFieldName
            diffTitleLabel = 'Main - Reference'

            remapObservationsSubtask = None

        for season in seasons:
            if refConfig is None:
                obsFileName = build_config_full_path(
                        config, 'seaIceObservations',
                        'thickness{}_{}'.format(hemisphere, season))

                remapObservationsSubtask = RemapObservedThickClimatology(
                        parentTask=self, seasons=[season],
                        fileName=obsFileName,
                        outFilePrefix='{}{}_{}'.format(refFieldName,
                                                       hemisphere,
                                                       season),
                        comparisonGridNames=comparisonGridNames,
                        subtaskName='remapObservations{}'.format(season))
                self.add_subtask(remapObservationsSubtask)

            for comparisonGridName in comparisonGridNames:

                imageDescription = \
                    '{} Climatology Map of {}-Hemisphere Sea-Ice ' \
                    'Thickness.'.format(season, hemisphereLong)
                imageCaption = imageDescription
                galleryGroup = \
                    '{}-Hemisphere Sea-Ice Thickness'.format(
                            hemisphereLong)
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                        self, hemisphere, season, comparisonGridName,
                        remapClimatologySubtask, remapObservationsSubtask,
                        refConfig)

                subtask.set_plot_info(
                        outFileLabel='icethick{}'.format(hemisphere),
                        fieldNameInTitle='Sea ice thickness',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'm',
                        imageDescription=imageDescription,
                        imageCaption=imageCaption,
                        galleryGroup=galleryGroup,
                        groupSubtitle=None,
                        groupLink='{}_thick'.format(hemisphere.lower()),
                        galleryName=galleryName,
                        maskValue=0)

                self.add_subtask(subtask)

        # }}}

    # }}}


class RemapObservedThickClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping sea ice thickness observations

    Authors
    -------
    Xylar Asay-Davis
    """

    def get_observation_descriptor(self, fileName):  # {{{
        '''
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid

        Authors
        -------
        Xylar Asay-Davis
        '''

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')
        return obsDescriptor  # }}}

    def build_observational_dataset(self, fileName):  # {{{
        '''
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

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_dataset(fileName)
        dsObs.rename({'HI': 'seaIceThick'}, inplace=True)
        return dsObs
        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
