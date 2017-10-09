
from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import os.path

import numpy.ma as ma
import numpy as np

import xarray as xr

from ..shared import AnalysisTask

from ..shared.climatology import get_comparison_descriptor, \
    get_remapper, get_observation_climatology_file_names, \
    remap_and_write_climatology, RemapMpasClimatologySubtask

from ..shared.grid import LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path
from ..shared.io import write_netcdf
from ..shared.html import write_image_xml


class ClimatologyMapSeaIce(AnalysisTask):  # {{{
    """
    An analysis task for comparison of 2D model fields against observations.

    Attributes
    ----------

    fieldName : str
        A short name of the field being analyzed

    fieldNameInTitle : str
        An equivalent name of the field appropriate for figure titles

    mpasFieldName : str
        The name of the MPAS timeSeriesStatsMonthly variable to be analyzed

    iselValues : dict
        A dictionary of dimensions and indices (or ``None``) used to extract
        a slice of the MPAS field.

    unitsLabel : str
        The units for the field being plotted.

    galleryGroup : str
        In the generated website, the name of the group of analysis galleries
        to which this analysis belongs

    groupLink : str
        the link within the component webpage used to identify the gallery
        group (a short version of the gallery group name with no spaces or
        punctuation)

    groupSubtitle : str
        a subtitle for the gallery group

    obsAndPlotInfo : list of dict
        A list of dictionaries containing information about each set of
        observations and each plot to be made.  Keys include ``season``,
        ``observationPrefix``, ``obsFileName``, ``observationTitleLabel``,
        ``outFileLabel``, ``galleryName``, ``imageDescription``,
        ``imageCaption``.

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    remapClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask that remaps the model climatologies this task will plot.
        The ``remapClimatologySubtask`` is needed to determine the file names
        of the climatology output.

    remapObservationsSubtask : ``RemapObservationsSubtask``
        The subtask that remaps the observational climatologies this task will
        plot.

    Authors
    -------
    Xylar Asay-Davis
    """

    def __init__(self, config, mpasClimatologyTask, taskName, tags):  # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        taskName :  str
            The name of the task, typically the same as the subclass name
            except starting with lowercase (e.g. 'climatologyMapSST' for class
            'ClimatologyMapSST')

        tags :  list of str
            Tags used to describe the task (e.g. 'climatology',
            'horizontalMap').  These are used to determine which tasks are
            generated (e.g. 'all_transect' or 'no_climatology' in the
            'generate' flags)

        Authors
        -------
        Xylar Asay-Davis

        '''
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIce, self).__init__(config=config,
                                                   taskName=taskName,
                                                   componentName='seaIce',
                                                   tags=tags)

        self.mpasClimatologyTask = mpasClimatologyTask

        sectionName = self.taskName

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        comparisonGridNames = ['latlon']

        if len(comparisonGridNames) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of comparison grids'.format(sectionName))

        # the variable self.mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        self.remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=self.mpasClimatologyTask,
            parentTask=self,
            climatologyName='{}{}'.format(self.fieldName, self.hemisphere),
            variableList=[self.mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=self.iselValues)

        # add a subtask for remapping the observations
        self.remapObservationsSubtask = RemapObservationsSubtask(
                self, seasons)
        self.add_subtask(self.remapObservationsSubtask)

        for info in self.obsAndPlotInfo:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(self, info)
                self.add_subtask(subtask)

        # }}}


class RemapObservationsSubtask(AnalysisTask):  # {{{
    """
    An analysis task for comparison of 2D model fields against observations.

    Attributes
    ----------
    parentTask :  ``ClimatologyMapSeaIce``
        The parent (master) task for this subtask

    seasons : list of str
        A list of seasons (keys in ``shared.constants.monthDictionary``)
        over which the climatology should be computed or ['none'] if only
        monthly climatologies are needed.

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, parentTask, seasons):  # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``ClimatologyMapOcean``
            The parent (master) task for this subtask

        seasons : list of str, optional
            A list of seasons (keys in ``shared.constants.monthDictionary``)
            to be computed or ['none'] (not ``None``) if only monthly
            climatologies are needed.

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.parentTask = parentTask
        self.seasons = seasons
        subtaskName = 'remapObservations'
        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(RemapObservationsSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='ocean', tags=tags)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(RemapObservationsSubtask, self).setup_and_check()

        # we set up the remapper here because ESFM_RegridWeightGen seems to
        # have trouble if it runs in another process (or in several at once)
        self._setup_obs_remapper()

        # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        config = self.config
        parentTask = self.parentTask

        hemisphere = parentTask.hemisphere

        for info in parentTask.obsAndPlotInfo:
            season = info['season']

            obsFileName = info['obsFileName']

            if not os.path.isfile(obsFileName):
                raise OSError('Obs file {} not found.'.format(
                    obsFileName))

            fieldName = '{}{}'.format(parentTask.obsFieldName, hemisphere)

            (obsClimatologyFileName, obsRemappedFileName) = \
                get_observation_climatology_file_names(
                    config=config, fieldName=fieldName,
                    monthNames=season, componentName=self.componentName,
                    remapper=self.obsRemapper)

            if not os.path.exists(obsRemappedFileName):

                # load the observations the first time
                seasonalClimatology = parentTask._build_observational_dataset(
                        obsFileName)
                write_netcdf(seasonalClimatology, obsClimatologyFileName)

                if self.obsRemapper.mappingFileName is None:
                    # no need to remap because the observations are on the
                    # comparison grid already
                    os.symlink(obsClimatologyFileName, obsRemappedFileName)
                else:
                    remap_and_write_climatology(
                        config, seasonalClimatology,
                        obsClimatologyFileName,
                        obsRemappedFileName, self.obsRemapper,
                        logger=self.logger)

        # }}}

    def _setup_obs_remapper(self):  # {{{
        """
        Set up the remapper for remapping from the MPAS to the comparison
        grid.

        Authors
        -------
        Xylar Asay-Davis
        """
        config = self.config
        parentTask = self.parentTask
        hemisphere = parentTask.hemisphere

        comparisonDescriptor = get_comparison_descriptor(
                config, comparisonGridName='latlon')

        info = parentTask.obsAndPlotInfo[0]

        obsFileName = info['obsFileName']

        obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')

        fieldName = '{}{}'.format(parentTask.obsFieldName, hemisphere)
        self.obsRemapper = get_remapper(
                config=config,
                sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(fieldName),
                method=config.get('seaIceObservations',
                                  'interpolationMethod'),
                logger=self.logger)
        # }}}
    # }}}


class PlotClimatologyMapSubtask(AnalysisTask):  # {{{
    """
    An analysis task for plotting 2D model fields against observations.

    Attributes
    ----------
    parentTask :  ``ClimatologyMapOcean``
        The parent (master) task for this subtask

    plotInfo : dict
        A dictionary containing information specific to the plot to be
        made.  Keys include ``season``, ``observationPrefix``,
        ``obsFileName``, ``observationTitleLabel``, ``outFileLabel``,
        ``galleryName``, ``imageDescription``, ``imageCaption``.

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
    """

    def __init__(self, parentTask, plotInfo):  # {{{
        '''
        Construct one analysis subtask for each plot (i.e. each season and
        comparison grid) and a subtask for computing climatologies.

        Parameters
        ----------
        parentTask :  ``ClimatologyMapOcean``
            The parent (master) task for this subtask

        plotInfo : dict
            A dictionary containing information specific to the plot to be
            made.  Keys include ``season``, ``observationPrefix``,
            ``obsFileName``, ``observationTitleLabel``, ``outFileLabel``,
            ``galleryName``, ``imageDescription``, ``imageCaption``.

        Authors
        -------
        Xylar Asay-Davis

        '''

        self.parentTask = parentTask
        self.plotInfo = plotInfo
        subtaskName = 'plot{}'.format(plotInfo['season'])
        if 'observationPrefix' in plotInfo.keys():
            subtaskName = '{}_{}'.format(subtaskName,
                                         plotInfo['observationPrefix'])
        config = parentTask.config
        taskName = parentTask.taskName
        tags = parentTask.tags

        # call the constructor from the base class (AnalysisTask)
        super(PlotClimatologyMapSubtask, self).__init__(
                config=config, taskName=taskName, subtaskName=subtaskName,
                componentName='ocean', tags=tags)

        # this task should not run until the remapping subtasks are done, since
        # it relies on data from those subtasks
        self.run_after(parentTask.remapClimatologySubtask)
        self.run_after(parentTask.remapObservationsSubtask)
        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup
        super(PlotClimatologyMapSubtask, self).setup_and_check()

        config = self.config
        self.startYear = config.getint('climatology', 'startYear')
        self.endYear = config.getint('climatology', 'endYear')
        self.startDate = config.get('climatology', 'startDate')
        self.endDate = config.get('climatology', 'endDate')

        mainRunName = config.get('runs', 'mainRunName')

        self.xmlFileNames = []

        info = self.plotInfo
        season = info['season']
        outFileLabel = info['outFileLabel']

        filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                outFileLabel, mainRunName,
                season, self.startYear, self.endYear)
        self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                    filePrefix))
        info['filePrefix'] = filePrefix

        # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        parentTask = self.parentTask
        info = self.plotInfo
        config = self.config
        season = info['season']
        fieldName = parentTask.fieldName
        if 'observationPrefix' in info.keys():
            self.logger.info("\nPlotting 2-d maps of {} climatologies for "
                             "{} from {} observations...".format(
                                     parentTask.fieldNameInTitle,
                                     season, info['observationPrefix']))
        else:
            self.logger.info("\nPlotting 2-d maps of {} climatologies for "
                             "{}...".format(parentTask.fieldNameInTitle,
                                            season))

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')
        startYear = self.startYear
        endYear = self.endYear

        comparisonDescriptor = get_comparison_descriptor(
                config, comparisonGridName='latlon')

        hemisphere = parentTask.hemisphere
        sectionName = parentTask.sectionName
        vertical = config.getboolean(sectionName, 'vertical')
        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        season = info['season']

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, sectionName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, sectionName, suffix='Difference')

        referenceLongitude = config.getfloat(sectionName,
                                             'referenceLongitude')
        minimumLatitude = config.getfloat(sectionName,
                                          'minimumLatitude')

        remappedFileName = parentTask.remapClimatologySubtask.get_file_name(
            season=season, stage='remapped', comparisonGridName='latlon')
        remappedClimatology = xr.open_dataset(remappedFileName)

        modelOutput = remappedClimatology[parentTask.mpasFieldName].values
        if parentTask.maskValue is not None:
            modelOutput = ma.masked_values(modelOutput, parentTask.maskValue)
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        obsFileName = info['obsFileName']

        obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')

        fieldName = '{}{}'.format(parentTask.obsFieldName, hemisphere)
        obsRemapper = get_remapper(
                config=config,
                sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(fieldName),
                method=config.get('seaIceObservations',
                                  'interpolationMethod'),
                logger=self.logger)

        (obsClimatologyFileName, obsRemappedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=fieldName,
                monthNames=season, componentName=self.componentName,
                remapper=obsRemapper)

        remappedClimatology = xr.open_dataset(obsRemappedFileName)

        observations = remappedClimatology[parentTask.obsFieldName].values
        if parentTask.maskValue is not None:
            observations = ma.masked_values(observations, parentTask.maskValue)

        difference = modelOutput - observations

        startYear = self.startYear
        endYear = self.endYear
        observationTitleLabel = info['observationTitleLabel']
        filePrefix = info['filePrefix']
        title = '{} ({}, years {:04d}-{:04d})'.format(
            parentTask.fieldNameInTitle, season, startYear, endYear)
        fileout = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
        plot_polar_comparison(
            config,
            lonTarg,
            latTarg,
            modelOutput,
            observations,
            difference,
            colormapResult,
            colorbarLevelsResult,
            colormapDifference,
            colorbarLevelsDifference,
            title=title,
            fileout=fileout,
            plotProjection=plotProjection,
            latmin=minimumLatitude,
            lon0=referenceLongitude,
            modelTitle=mainRunName,
            obsTitle=observationTitleLabel,
            diffTitle='Model-Observations',
            cbarlabel=parentTask.unitsLabel,
            vertical=vertical)

        galleryName = info['galleryName']
        imageDescription = info['imageDescription']
        imageCaption = info['imageCaption']
        write_image_xml(
            config,
            filePrefix,
            componentName='Sea Ice',
            componentSubdirectory='sea_ice',
            galleryGroup=parentTask.galleryGroup,
            groupSubtitle=parentTask.groupSubtitle,
            groupLink=parentTask.groupLink,
            gallery=galleryName,
            thumbnailDescription=season,
            imageDescription=imageDescription,
            imageCaption=imageCaption)
        # }}}

    # }}}


class ClimatologyMapSeaIceConc(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice concentration against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, mpasClimatologyTask, hemisphere):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        Authors
        -------
        Xylar Asay-Davis
        """
        taskName = 'climatologyMapSeaIceConc{}'.format(hemisphere)

        self.sectionName = taskName

        self.fieldName = 'seaIceConc'
        self.hemisphere = hemisphere
        self.fieldNameInTitle = 'Sea ice concentration'
        self.seasons = config.getExpression(self.sectionName, 'seasons')
        self.observationPrefixes = config.getExpression(self.sectionName,
                                                        'observationPrefixes')

        self.mpasFieldName = 'timeMonthly_avg_iceAreaCell'
        self.iselValues = None

        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        self.obsAndPlotInfo = []
        for prefix in self.observationPrefixes:
            for season in self.seasons:
                localDict = {}
                localDict['season'] = season
                localDict['observationPrefix'] = prefix
                localDict['obsFileName'] = \
                    build_config_full_path(
                        config, 'seaIceObservations',
                        'concentration{}{}_{}'.format(prefix,
                                                      hemisphere,
                                                      season))
                localDict['observationTitleLabel'] = \
                    'Observations (SSM/I {})'.format(prefix)
                localDict['outFileLabel'] = \
                    'iceconc{}{}'.format(prefix, hemisphere)

                localDict['galleryName'] = \
                    'Observations: SSM/I {}'.format(prefix)

                localDict['imageDescription'] = \
                    '{} Climatology Map of {}-Hemisphere Sea-Ice ' \
                    'Concentration'.format(season, hemisphereLong)

                localDict['imageCaption'] = \
                    '{}. <br> Observations: SSM/I {}'.format(
                            localDict['imageDescription'], prefix)

                self.obsAndPlotInfo.append(localDict)

        self.obsFieldName = 'AICE'

        self.unitsLabel = 'fraction'

        self.maskValue = None

        # variables for XML and webpages
        self.galleryGroup = '{}-Hemisphere Sea-Ice Concentration'.format(
                 hemisphereLong)
        self.groupSubtitle = None
        self.groupLink = '{}_conc'.format(hemisphere.lower())

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(
                config=config, mpasClimatologyTask=mpasClimatologyTask,
                taskName=taskName, tags=tags)

        # }}}

    def _build_observational_dataset(self, obsFileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(obsFileName)

        return dsObs  # }}}

    # }}}


class ClimatologyMapSeaIceThick(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice thickness against observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, mpasClimatologyTask, hemisphere):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        hemisphere : {'NH', 'SH'}
            The hemisphere to plot

        Authors
        -------
        Xylar Asay-Davis
        """
        taskName = 'climatologyMapSeaIceThick{}'.format(hemisphere)

        self.sectionName = taskName

        self.fieldName = 'seaIceThick'
        self.hemisphere = hemisphere
        self.fieldNameInTitle = 'Sea ice thickness'
        self.seasons = config.getExpression(self.sectionName, 'seasons')

        self.mpasFieldName = 'timeMonthly_avg_iceVolumeCell'
        self.iselValues = None

        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        self.obsAndPlotInfo = []
        for season in self.seasons:
            localDict = {}
            localDict['season'] = season
            localDict['obsFileName'] = \
                build_config_full_path(
                    config, 'seaIceObservations',
                    'thickness{}_{}'.format(hemisphere, season))
            localDict['observationTitleLabel'] = 'Observations (ICESat)'
            localDict['outFileLabel'] = 'icethick{}'.format(hemisphere)

            localDict['galleryName'] = 'Observations: ICESat'

            localDict['imageDescription'] = \
                '{} Climatology Map of {}-Hemisphere Sea-Ice ' \
                'Thickness.'.format(season, hemisphereLong)

            localDict['imageCaption'] = localDict['imageDescription']

            self.obsAndPlotInfo.append(localDict)

        self.obsFieldName = 'HI'

        self.unitsLabel = 'm'

        self.maskValue = 0

        # variables for XML and webpages
        self.galleryGroup = '{}-Hemisphere Sea-Ice Thicknesss'.format(
                 hemisphereLong)
        self.groupSubtitle = None
        self.groupLink = '{}_thick'.format(hemisphere.lower())

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(
                config=config, mpasClimatologyTask=mpasClimatologyTask,
                taskName=taskName, tags=tags)

        # }}}

    def _build_observational_dataset(self, obsFileName):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(obsFileName)

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
