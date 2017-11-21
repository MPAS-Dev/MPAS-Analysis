import os
import os.path

import numpy.ma as ma
import numpy as np

import xarray as xr

from ..shared.climatology import get_comparison_descriptor, \
    get_remapper, get_observation_climatology_file_names, \
    remap_and_write_climatology, RemapMpasClimatologySubtask

from ..shared.grid import LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path, make_directories
from ..shared.io import write_netcdf
from ..shared.html import write_image_xml

from .sea_ice_analysis_task import SeaIceAnalysisTask


class ClimatologyMapSeaIce(SeaIceAnalysisTask):
    """
    General comparison of 2-d model fields against data.  Currently only
    supports sea ice concentration and sea ice thickness

    Attributes
    ----------
    remapClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask that remaps the climatologies this task will plot.
        The ``remapClimatologySubtask`` is needed to determine the file names
        of the climatology output.

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted


    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
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

        # the variable self.mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        self.remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=self.mpasClimatologyTask,
            parentTask=self,
            climatologyName=self.fieldName,
            variableList=[self.mpasFieldName],
            seasons=seasons,
            iselValues=self.iselValues)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (SeaIceAnalysisTask),
        # which will perform some common setup
        super(ClimatologyMapSeaIce, self).setup_and_check()

        self.startYear = self.mpasClimatologyTask.startYear
        self.startDate = self.mpasClimatologyTask.startDate
        self.endYear = self.mpasClimatologyTask.endYear
        self.endDate = self.mpasClimatologyTask.endDate

        mainRunName = self.config.get('runs', 'mainRunName')

        # we set up the remapper here because ESFM_RegridWeightGen seems to
        # have trouble if it runs in another process (or in several at once)
        self._setup_obs_remapper()

        self.xmlFileNames = []

        for info in self.obsAndPlotInfo:
            season = info['season']
            outFileLabel = info['outFileLabel']

            filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                    outFileLabel, mainRunName,
                    season, self.startYear, self.endYear)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            info['filePrefix'] = filePrefix

        # make the mapping directory, because doing so within each process
        # seems to be giving ESMF_RegridWeightGen some trouble
        mappingSubdirectory = \
            build_config_full_path(self.config, 'output',
                                   'mappingSubdirectory')
        make_directories(mappingSubdirectory)

        # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        self.logger.info("\nPlotting 2-d maps of {} climatologies...".format(
            self.fieldNameInTitle))

        config = self.config

        mainRunName = config.get('runs', 'mainRunName')
        startYear = self.startYear
        endYear = self.endYear

        hemisphere = self.hemisphere
        sectionName = self.sectionName
        vertical = config.getboolean(sectionName, 'vertical')
        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        for info in self.obsAndPlotInfo:
            season = info['season']

            (colormapResult, colorbarLevelsResult) = setup_colormap(
                config, sectionName, suffix='Result')
            (colormapDifference, colorbarLevelsDifference) = setup_colormap(
                config, sectionName, suffix='Difference')

            referenceLongitude = config.getfloat(sectionName,
                                                 'referenceLongitude')
            minimumLatitude = config.getfloat(sectionName,
                                              'minimumLatitude')

            fieldName = '{}{}'.format(self.mpasFieldName, hemisphere)

            remappedFileName = self.remapClimatologySubtask.get_file_name(
                season=season, stage='remapped', comparisonGridName='latlon')
            remappedClimatology = xr.open_dataset(remappedFileName)

            modelOutput = remappedClimatology[self.mpasFieldName].values
            if self.maskValue is not None:
                modelOutput = ma.masked_values(modelOutput, self.maskValue)
            lon = remappedClimatology['lon'].values
            lat = remappedClimatology['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            obsFileName = info['obsFileName']

            if not os.path.isfile(obsFileName):
                raise OSError('Obs file {} not found.'.format(
                    obsFileName))

            fieldName = '{}{}'.format(self.obsFieldName, hemisphere)

            (obsClimatologyFileName, obsRemappedFileName) = \
                get_observation_climatology_file_names(
                    config=config, fieldName=fieldName,
                    monthNames=season, componentName=self.componentName,
                    remapper=self.obsRemapper)

            if not os.path.exists(obsRemappedFileName):

                # load the observations the first time
                seasonalClimatology = self._build_observational_dataset(
                        obsFileName)
                write_netcdf(seasonalClimatology, obsClimatologyFileName)

                if self.obsRemapper is None:
                    # no need to remap because the observations are on the
                    # comparison grid already
                    remappedClimatology = seasonalClimatology
                else:
                    remappedClimatology = \
                        remap_and_write_climatology(
                            config, seasonalClimatology,
                            obsClimatologyFileName,
                            obsRemappedFileName, self.obsRemapper,
                            logger=self.logger)

            else:

                remappedClimatology = xr.open_dataset(obsRemappedFileName)

            observations = remappedClimatology[self.obsFieldName].values
            if self.maskValue is not None:
                observations = ma.masked_values(observations, self.maskValue)

            difference = modelOutput - observations

            startYear = self.startYear
            endYear = self.endYear
            observationTitleLabel = info['observationTitleLabel']
            filePrefix = info['filePrefix']
            title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, season, startYear, endYear)
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
                cbarlabel=self.unitsLabel,
                vertical=vertical)

            galleryName = info['galleryName']
            imageDescription = info['imageDescription']
            imageCaption = info['imageCaption']
            write_image_xml(
                config,
                filePrefix,
                componentName='Sea Ice',
                componentSubdirectory='sea_ice',
                galleryGroup=self.galleryGroup,
                groupSubtitle=self.groupSubtitle,
                groupLink=self.groupLink,
                gallery=galleryName,
                thumbnailDescription=season,
                imageDescription=imageDescription,
                imageCaption=imageCaption)
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
        hemisphere = self.hemisphere

        info = self.obsAndPlotInfo[0]

        fieldName = '{}{}'.format(self.mpasFieldName, hemisphere)

        obsFileName = info['obsFileName']

        obsDescriptor = LatLonGridDescriptor.read(fileName=obsFileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')

        comparisonDescriptor = get_comparison_descriptor(
                config=config, comparisonGridName='latlon')

        fieldName = '{}{}'.format(self.obsFieldName, hemisphere)
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

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(
                config=config, mpasClimatologyTask=mpasClimatologyTask,
                taskName=taskName, tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        hemisphere = self.hemisphere
        if hemisphere == 'NH':
            hemisphereLong = 'Northern'
        else:
            hemisphereLong = 'Southern'

        self.obsAndPlotInfo = []
        for prefix in self.observationPrefixes:
            for season in self.seasons:
                localDict = {}
                localDict['season'] = season
                localDict['obsFileName'] = \
                    build_config_full_path(
                        self.config, 'seaIceObservations',
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

        # call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceConc, self).setup_and_check()

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

        tags = ['climatology', 'horizontalMap', self.fieldName]

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(
                config=config, mpasClimatologyTask=mpasClimatologyTask,
                taskName=taskName, tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        hemisphere = self.hemisphere

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
                    self.config, 'seaIceObservations',
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

        # call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceThick, self).setup_and_check()

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
