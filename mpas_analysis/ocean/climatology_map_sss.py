from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr
import datetime

from ..shared import AnalysisTask

from ..shared.io.utility import build_config_full_path

from ..shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask, RemapMpasReferenceClimatologySubtask

from .plot_climatology_map_subtask import PlotClimatologyMapSubtask

from ..shared.grid import LatLonGridDescriptor


class ClimatologyMapSSS(AnalysisTask):  # {{{
    """
    An analysis task for comparison of sea surface salinity (sss) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, mpasClimatologyTask,
                 mpasRefClimatologyTask=None):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        mpasRefClimatologyTask : ``MpasClimatologyTask``, optional
            The task that produced the climatology from a reference run to be
            remapped and plotted, including anomalies with respect to the main
            run

        Authors
        -------
        Xylar Asay-Davis
        """
        fieldName = 'sss'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSSS, self).__init__(
                config=config, taskName='climatologyMapSSS',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', fieldName])

        mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        iselValues = {'nVertLevels': 0}

        sectionName = self.taskName

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
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if mpasRefClimatologyTask is None:

            refTitleLabel = \
                'Observations (Aquarius, 2011-2014)'

            observationsDirectory = build_config_full_path(
                config, 'oceanObservations',
                '{}Subdirectory'.format(fieldName))

            obsFileName = \
                "{}/Aquarius_V3_SSS_Monthly.nc".format(
                    observationsDirectory)
            refFieldName = 'sss'
            outFileLabel = 'sssAquarius'
            galleryName = 'Observations: Aquarius'

            remapObservationsSubtask = RemapObservedSSSClimatology(
                    parentTask=self, seasons=seasons, fileName=obsFileName,
                    outFilePrefix=refFieldName,
                    comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            remapRefClimatologySubtask = None
            diffTitleLabel = 'Model - Observations'

        else:
            remapRefClimatologySubtask = RemapMpasReferenceClimatologySubtask(
                mpasClimatologyTask=mpasRefClimatologyTask,
                parentTask=self,
                climatologyName=fieldName,
                variableList=[mpasFieldName],
                comparisonGridNames=comparisonGridNames,
                seasons=seasons,
                iselValues=iselValues)
            remapObservationsSubtask = None
            refRunName = mpasRefClimatologyTask.config.get(
                    'runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Ref: {}'.format(refRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'sss'
            diffTitleLabel = 'Main - Reference'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(self, season,
                                                    comparisonGridName,
                                                    remapClimatologySubtask,
                                                    remapObservationsSubtask,
                                                    remapRefClimatologySubtask)

                subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='SSS',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'PSU',
                        imageCaption='Mean Sea Surface Salinity',
                        galleryGroup='Sea Surface Salinity',
                        groupSubtitle=None,
                        groupLink='sss',
                        galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}
    # }}}


class RemapObservedSSSClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping SSS observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
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
                                                  latVarName='lat',
                                                  lonVarName='lon')
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

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        dsObs = xr.open_dataset(fileName)
        dsObs.rename({'time': 'Time', 'SSS': 'sss'}, inplace=True)
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
