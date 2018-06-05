# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.sea_ice.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.grid import LatLonGridDescriptor


class ClimatologyMapSeaIceConc(AnalysisTask):  # {{{
    """
    An analysis task for comparison of sea ice concentration against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

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
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        taskName = 'climatologyMapSeaIceConc{}'.format(hemisphere)

        fieldName = 'seaIceConc'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(
                config=config, taskName=taskName,
                componentName='seaIce',
                tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

        mpasFieldName = 'timeMonthly_avg_iceAreaCell'
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
            self._add_obs_tasks(seasons, comparisonGridNames, hemisphere,
                                hemisphereLong, remapClimatologySubtask,
                                mpasFieldName)
        else:
            self._add_ref_tasks(seasons, comparisonGridNames, hemisphere,
                                hemisphereLong, remapClimatologySubtask,
                                refConfig, mpasFieldName,
                                fieldName, iselValues)
        # }}}

    def _add_obs_tasks(self, seasons, comparisonGridNames, hemisphere,
                       hemisphereLong, remapClimatologySubtask,
                       mpasFieldName):  # {{{
        config = self.config
        obsFieldName = 'seaIceConc'
        sectionName = self.taskName

        observationPrefixes = config.getExpression(sectionName,
                                                   'observationPrefixes')
        for prefix in observationPrefixes:
            for season in seasons:
                observationTitleLabel = \
                    'Observations (SSM/I {})'.format(prefix)

                obsFileName = build_config_full_path(
                        config=config, section='seaIceObservations',
                        relativePathOption='concentration{}{}_{}'.format(
                                prefix, hemisphere, season),
                        relativePathSection=sectionName)

                remapObservationsSubtask = RemapObservedConcClimatology(
                        parentTask=self, seasons=[season],
                        fileName=obsFileName,
                        outFilePrefix='{}{}{}_{}'.format(
                                obsFieldName,  prefix, hemisphere, season),
                        comparisonGridNames=comparisonGridNames,
                        subtaskName='remapObservations_{}{}'.format(
                                prefix,  season))
                self.add_subtask(remapObservationsSubtask)
                for comparisonGridName in comparisonGridNames:

                    imageDescription = \
                        '{} Climatology Map of {}-Hemisphere Sea-Ice ' \
                        'Concentration'.format(season, hemisphereLong)
                    imageCaption = \
                        '{}. <br> Observations: SSM/I {}'.format(
                            imageDescription, prefix)
                    galleryGroup = \
                        '{}-Hemisphere Sea-Ice Concentration'.format(
                                hemisphereLong)
                    # make a new subtask for this season and comparison
                    # grid
                    subtask = PlotClimatologyMapSubtask(
                        self, hemisphere, season, comparisonGridName,
                        remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        subtaskSuffix=prefix)

                    subtask.set_plot_info(
                        outFileLabel='iceconc{}{}'.format(prefix,
                                                          hemisphere),
                        fieldNameInTitle='Sea ice concentration',
                        mpasFieldName=mpasFieldName,
                        refFieldName=obsFieldName,
                        refTitleLabel=observationTitleLabel,
                        diffTitleLabel='Model - Observations',
                        unitsLabel=r'fraction',
                        imageDescription=imageDescription,
                        imageCaption=imageCaption,
                        galleryGroup=galleryGroup,
                        groupSubtitle=None,
                        groupLink='{}_conc'.format(hemisphere.lower()),
                        galleryName='Observations: SSM/I {}'.format(
                                prefix))

                    self.add_subtask(subtask)
        # }}}

    def _add_ref_tasks(self, seasons, comparisonGridNames, hemisphere,
                       hemisphereLong, remapClimatologySubtask,
                       refConfig, mpasFieldName, fieldName,
                       iselValues):  # {{{

        refRunName = refConfig.get('runs', 'mainRunName')
        galleryName = None
        refTitleLabel = 'Ref: {}'.format(refRunName)

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
                    self, hemisphere, season, comparisonGridName,
                    remapClimatologySubtask, refConfig=refConfig)

                subtask.set_plot_info(
                        outFileLabel='iceconc{}'.format(hemisphere),
                        fieldNameInTitle='Sea ice concentration',
                        mpasFieldName=mpasFieldName,
                        refFieldName=mpasFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel='Main - Reference',
                        unitsLabel=r'fraction',
                        imageDescription=imageDescription,
                        imageCaption=imageCaption,
                        galleryGroup=galleryGroup,
                        groupSubtitle=None,
                        groupLink='{}_conc'.format(hemisphere.lower()),
                        galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}
    # }}}


class RemapObservedConcClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping sea ice concentration observations
    """
    # Authors
    # -------
    # Xylar Asay-Davis

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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        dsObs.rename({'AICE': 'seaIceConc'}, inplace=True)
        return dsObs
        # }}}
    # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
