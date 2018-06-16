from __future__ import absolute_import, division, print_function, \
    unicode_literals

import xarray as xr

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask

from mpas_analysis.shared.grid import LatLonGridDescriptor

from mpas_analysis.shared.constants import constants


class ClimatologyMapSSH(AnalysisTask):  # {{{
    """
    An analysis task for comparison of sea surface height (ssh) against
    observations
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

        fieldName = 'ssh'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSSH, self).__init__(
                config=config, taskName='climatologyMapSSH',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

        mpasFieldName = 'timeMonthly_avg_pressureAdjustedSSH'

        iselValues = None

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

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        remapClimatologySubtask = RemapSSHClimatology(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if refConfig is None:

            refTitleLabel = 'Observations (AVISO Dynamic ' \
                'Topography, 1993-2010)'

            observationsDirectory = build_config_full_path(
                config, 'oceanObservations',
                '{}Subdirectory'.format(fieldName))

            obsFileName = \
                "{}/zos_AVISO_L4_199210-201012.nc".format(
                    observationsDirectory)
            refFieldName = 'zos'
            outFileLabel = 'sshAVISO'
            galleryName = 'Observations: AVISO'

            remapObservationsSubtask = RemapObservedSSHClimatology(
                    parentTask=self, seasons=seasons, fileName=obsFileName,
                    outFilePrefix=refFieldName,
                    comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            refRunName = refConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Ref: {}'.format(refRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'ssh'
            diffTitleLabel = 'Main - Reference'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(self, season,
                                                    comparisonGridName,
                                                    remapClimatologySubtask,
                                                    remapObservationsSubtask,
                                                    refConfig, removeMean=True)

                subtask.set_plot_info(
                        outFileLabel=outFileLabel,
                        fieldNameInTitle='Zero-mean SSH',
                        mpasFieldName=mpasFieldName,
                        refFieldName=refFieldName,
                        refTitleLabel=refTitleLabel,
                        diffTitleLabel=diffTitleLabel,
                        unitsLabel=r'cm',
                        imageCaption='Mean Sea Surface Height',
                        galleryGroup='Sea Surface Height',
                        groupSubtitle=None,
                        groupLink='ssh',
                        galleryName=galleryName)

                self.add_subtask(subtask)
        # }}}
    # }}}


class RemapSSHClimatology(RemapMpasClimatologySubtask):  # {{{
    """
    Change units from m to cm
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def customize_masked_climatology(self, climatology, season):  # {{{
        """
        Mask the melt rates using ``landIceMask`` and rescale it to m/yr

        Parameters
        ----------
        climatology : ``xarray.Dataset```
            The MPAS climatology data set that has had a mask added but has
            not yet been remapped

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

        fieldName = self.variableList[0]

        # scale the field to cm from m
        climatology[fieldName] = constants.cm_per_m * climatology[fieldName]

        return climatology  # }}}

    # }}}


class RemapObservedSSHClimatology(RemapObservedClimatologySubtask):  # {{{
    """
    A subtask for reading and remapping SSH observations
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
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        dsObs = xr.open_dataset(fileName)
        dsObs.rename({'time': 'Time'}, inplace=True)
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        # scale the field to cm from m
        dsObs['zos'] = constants.cm_per_m * dsObs['zos']

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
