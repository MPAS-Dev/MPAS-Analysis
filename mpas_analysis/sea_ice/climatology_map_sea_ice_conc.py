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

import xarray as xr
from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.shared.plot import PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_obs_path


class ClimatologyMapSeaIceConc(AnalysisTask):
    """
    An analysis task for comparison of sea ice concentration against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasClimatologyTask, hemisphere,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        taskName = 'climatologyMapSeaIceConc{}'.format(hemisphere)

        fieldName = 'seaIceConc'

        tags = ['climatology', 'horizontalMap', fieldName, 'publicObs']
        if hemisphere == 'NH':
            tags = tags + ['arctic']
        else:
            tags = tags + ['antarctic']

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(
            config=config, taskName=taskName,
            componentName='seaIce',
            tags=tags)

        mpasFieldName = 'timeMonthly_avg_iceAreaCell'
        iselValues = None

        sectionName = taskName

        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        # read in what seasons we want to plot
        seasons = config.getexpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = config.getexpression(sectionName,
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

        if controlConfig is None:
            self._add_obs_tasks(seasons, comparisonGridNames, hemisphere,
                                hemisphereLong, remapClimatologySubtask,
                                mpasFieldName)
        else:
            self._add_ref_tasks(seasons, comparisonGridNames, hemisphere,
                                hemisphereLong, remapClimatologySubtask,
                                controlConfig, mpasFieldName,
                                fieldName, iselValues)

    def _add_obs_tasks(self, seasons, comparisonGridNames, hemisphere,
                       hemisphereLong, remapClimatologySubtask,
                       mpasFieldName):
        config = self.config
        obsFieldName = 'seaIceConc'
        sectionName = self.taskName

        observationPrefixes = config.getexpression(sectionName,
                                                   'observationPrefixes')
        for prefix in observationPrefixes:
            for season in seasons:
                observationTitleLabel = \
                    'Observations (SSM/I {})'.format(prefix)

                obsFileName = build_obs_path(
                    config, 'seaIce',
                    relativePathOption='concentration{}{}_{}'.format(
                        prefix, hemisphere, season),
                    relativePathSection=sectionName)

                remapObservationsSubtask = RemapObservedConcClimatology(
                    parentTask=self, seasons=[season],
                    fileName=obsFileName,
                    outFilePrefix='{}{}{}_{}'.format(
                        obsFieldName, prefix, hemisphere, season),
                    comparisonGridNames=comparisonGridNames,
                    subtaskName='remapObservations_{}{}'.format(
                        prefix, season))
                self.add_subtask(remapObservationsSubtask)
                for comparisonGridName in comparisonGridNames:

                    imageDescription = \
                        'Climatology Map of {}-Hemisphere Sea-Ice ' \
                        'Concentration'.format(hemisphereLong)
                    imageCaption = \
                        '{}. <br> Observations: SSM/I {}'.format(
                            imageDescription, prefix)
                    galleryGroup = \
                        '{}-Hemisphere Sea-Ice Concentration'.format(
                            hemisphereLong)
                    # make a new subtask for this season and comparison
                    # grid

                    subtaskName = f'plot{season}_{comparisonGridName}_{prefix}'

                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self, season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        subtaskName=subtaskName)

                    subtask.set_plot_info(
                        outFileLabel='iceconc{}{}'.format(prefix,
                                                          hemisphere),
                        fieldNameInTitle='Sea ice concentration',
                        mpasFieldName=mpasFieldName,
                        refFieldName=obsFieldName,
                        refTitleLabel=observationTitleLabel,
                        diffTitleLabel='Model - Observations',
                        unitsLabel=r'fraction',
                        imageCaption=imageCaption,
                        galleryGroup=galleryGroup,
                        groupSubtitle=None,
                        groupLink='{}_conc'.format(hemisphere.lower()),
                        galleryName='Observations: SSM/I {}'.format(
                            prefix))

                    self.add_subtask(subtask)

    def _add_ref_tasks(self, seasons, comparisonGridNames, hemisphere,
                       hemisphereLong, remapClimatologySubtask,
                       controlConfig, mpasFieldName, fieldName,
                       iselValues):

        controlRunName = controlConfig.get('runs', 'mainRunName')
        galleryName = None
        refTitleLabel = 'Control: {}'.format(controlRunName)

        for season in seasons:
            for comparisonGridName in comparisonGridNames:

                imageDescription = \
                    '{} Climatology Map of {}-Hemisphere Sea-Ice ' \
                    'Concentration'.format(season, hemisphereLong)
                imageCaption = imageDescription
                galleryGroup = \
                    '{}-Hemisphere Sea-Ice Concentration'.format(
                        hemisphereLong)
                # make a new subtask for this season and comparison
                # grid
                subtask = PlotClimatologyMapSubtask(
                    parentTask=self, season=season,
                    comparisonGridName=comparisonGridName,
                    remapMpasClimatologySubtask=remapClimatologySubtask,
                    controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel='iceconc{}'.format(hemisphere),
                    fieldNameInTitle='Sea ice concentration',
                    mpasFieldName=mpasFieldName,
                    refFieldName=mpasFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel='Main - Control',
                    unitsLabel=r'fraction',
                    imageCaption=imageCaption,
                    galleryGroup=galleryGroup,
                    groupSubtitle=None,
                    groupLink='{}_conc'.format(hemisphere.lower()),
                    galleryName=galleryName)

                self.add_subtask(subtask)


class RemapObservedConcClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping sea ice concentration observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def get_observation_descriptor(self, fileName):
        """
        get a MeshDescriptor for the observation grid

        Parameters
        ----------
        fileName : str
            observation file name describing the source grid

        Returns
        -------
        obsDescriptor : ``MeshDescriptor``
            The descriptor for the observation grid
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # create a descriptor of the observation grid using the lat/lon
        # coordinates
        obsDescriptor = LatLonGridDescriptor.read(fileName=fileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')
        return obsDescriptor

    def build_observational_dataset(self, fileName):
        """
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
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'AICE': 'seaIceConc'})
        return dsObs
