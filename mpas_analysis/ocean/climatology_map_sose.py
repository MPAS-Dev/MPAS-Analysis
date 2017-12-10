'''
Analysis tasks for comparing Antarctic climatology maps against observations
and reanalysis data.

Authors
-------
Xylar Asay-Davis
'''

import xarray as xr

from ..shared import AnalysisTask

from .remap_depth_slices_subtask import RemapDepthSlicesSubtask

from .plot_climatology_map_subtask import PlotClimatologyMapSubtask

from ..shared.io.utility import build_config_full_path

from ..shared.climatology import RemapObservedClimatologySubtask, \
    get_antarctic_stereographic_projection

from ..shared.grid import ProjectionGridDescriptor

from ..shared.mpas_xarray import mpas_xarray


class ClimatologyMapSoseTemperature(AnalysisTask):  # {{{
    """
    An analysis task for comparison of antarctic temperature against SOSE
    fields

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, mpasClimatologyTask):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        Authors
        -------
        Xylar Asay-Davis
        """
        fieldName = 'temperature'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSoseTemperature, self).__init__(
                config=config, taskName='climatologyMapSoseTemperature',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', 'sose', fieldName])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        iselValues = None

        observationTitleLabel = 'State Estimate (SOSE)'

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'soseSubdirectory')

        obsFileName = \
            '{}/SOSE_2005-2010_monthly_pot_temp_6000.0x' \
            '6000.0km_10.0km_Antarctic_stereo.nc'.format(observationsDirectory)
        obsFieldName = 'theta'

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

        depths = config.getExpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of depths'.format(sectionName))

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        remapObservationsSubtask = RemapSoseClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=obsFieldName,
                fieldName=obsFieldName,
                botFieldName='botTheta',
                depths=depths,
                comparisonGridNames=comparisonGridNames)
        self.add_subtask(remapObservationsSubtask)
        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        depth=depth)

                    subtask.set_plot_info(
                        outFileLabel='tempSOSE',
                        fieldNameInTitle='Temperature',
                        mpasFieldName=mpasFieldName,
                        obsFieldName=obsFieldName,
                        observationTitleLabel=observationTitleLabel,
                        diffTitleLabel='Model - State Estimate',
                        unitsLabel=r'$^\circ$C',
                        imageCaption='Temperature',
                        galleryGroup='Temperature',
                        groupSubtitle=None,
                        groupLink='temp',
                        galleryName='State Estimate: SOSE')

                    self.add_subtask(subtask)
        # }}}

    # }}}


class ClimatologyMapSoseSalinity(AnalysisTask):  # {{{
    """
    An analysis task for comparison of antarctic salinity against SOSE
    fields

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, mpasClimatologyTask):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        Authors
        -------
        Xylar Asay-Davis
        """
        fieldName = 'salinity'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSoseSalinity, self).__init__(
                config=config, taskName='climatologyMapSoseSalinity',
                componentName='ocean',
                tags=['climatology', 'horizontalMap', 'sose', fieldName])

        sectionName = self.taskName

        mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        iselValues = None

        observationTitleLabel = 'State Estimate (SOSE)'

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', 'soseSubdirectory')

        obsFileName = \
            '{}/SOSE_2005-2010_monthly_salinity_6000.0x' \
            '6000.0km_10.0km_Antarctic_stereo.nc'.format(observationsDirectory)
        obsFieldName = 'salinity'

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

        depths = config.getExpression(sectionName, 'depths')

        if len(depths) == 0:
            raise ValueError('config section {} does not contain valid '
                             'list of depths'.format(sectionName))

        # the variable 'timeMonthly_avg_landIceFreshwaterFlux' will be added to
        # mpasClimatologyTask along with the seasons.
        remapClimatologySubtask = RemapDepthSlicesSubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            seasons=seasons,
            depths=depths,
            comparisonGridNames=comparisonGridNames,
            iselValues=iselValues)

        remapObservationsSubtask = RemapSoseClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=obsFieldName,
                fieldName=obsFieldName,
                botFieldName='botSalinity',
                depths=depths,
                comparisonGridNames=comparisonGridNames)
        self.add_subtask(remapObservationsSubtask)
        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                for depth in depths:
                    subtask = PlotClimatologyMapSubtask(
                        parentTask=self,
                        season=season,
                        comparisonGridName=comparisonGridName,
                        remapMpasClimatologySubtask=remapClimatologySubtask,
                        remapObsClimatologySubtask=remapObservationsSubtask,
                        depth=depth)

                    subtask.set_plot_info(
                        outFileLabel='salinSOSE',
                        fieldNameInTitle='Salinity',
                        mpasFieldName=mpasFieldName,
                        obsFieldName=obsFieldName,
                        observationTitleLabel=observationTitleLabel,
                        diffTitleLabel='Model - State Estimate',
                        unitsLabel=r'PSU',
                        imageCaption='Salinity',
                        galleryGroup='Salinity',
                        groupSubtitle=None,
                        groupLink='salin',
                        galleryName='State Estimate: SOSE')

                    self.add_subtask(subtask)
        # }}}

    # }}}


class RemapSoseClimatology(RemapObservedClimatologySubtask):
    # {{{
    """
    A subtask for reading and remapping SOSE fields to the comparison grid

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, parentTask, seasons, fileName, outFilePrefix,
                 fieldName, botFieldName, depths,
                 comparisonGridNames=['latlon'],
                 subtaskName='remapObservations'):
        # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``AnalysisTask``
            The parent (master) task for this subtask

        seasons : list of str
           A list of seasons (keys in ``constants.monthDictionary``) over
           which the climatology should be computed.

        fileName : str
            The name of the observation file

        outFilePrefix : str
            The prefix in front of output files and mapping files, typically
            the name of the field being remapped

        fieldName : str
            The name of the 3D field to remap

        botFieldName : str
            The name of the same field as ``fieldName`` but sampled at the
            sea floor

        depths : list of {None, float, 'top', 'bot'}
            A list of depths at which the climatology will be sliced in the
            vertical.

        comparisonGridNames : list of {'latlon', 'antarctic'}, optional
            The name(s) of the comparison grid to use for remapping.

        subtaskName : str, optional
            The name of the subtask

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.fieldName = fieldName
        self.botFieldName = botFieldName
        self.depths = depths

        # call the constructor from the base class (AnalysisTask)
        super(RemapSoseClimatology, self).__init__(
                parentTask, seasons, fileName, outFilePrefix,
                comparisonGridNames, subtaskName)
        # }}}

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

        # create a descriptor of the observation grid using the x/y polar
        # stereographic coordinates
        projection = get_antarctic_stereographic_projection()
        obsDescriptor = ProjectionGridDescriptor.read(
            projection, fileName=fileName, xVarName='x', yVarName='y')
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

        # Load MLD observational data
        dsObs = xr.open_dataset(fileName)

        dsObs = mpas_xarray.subset_variables(dsObs, [self.fieldName,
                                                     self.botFieldName,
                                                     'month', 'year'])
        slices = []
        field = dsObs[self.fieldName]
        botField = dsObs[self.botFieldName]
        for depth in self.depths:
            if depth == 'top':
                slices.append(field.sel(method='nearest', depth=0.).drop(
                        'depth'))
            elif depth == 'bot':
                slices.append(botField)
            else:
                slices.append(field.sel(method='nearest', depth=depth).drop(
                        'depth'))

        depthNames = [str(depth) for depth in self.depths]
        field = xr.concat(slices, dim='depthSlice')

        dsObs = xr.Dataset(data_vars={self.fieldName: field},
                           coords={'depthSlice': depthNames})

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
