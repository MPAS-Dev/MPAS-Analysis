import os
import os.path

import numpy.ma as ma
import numpy as np

import xarray as xr

from ..shared.climatology import get_remapper, \
    get_observation_climatology_file_names, \
    remap_and_write_climatology, \
    get_comparison_descriptor, MpasClimatology
from ..shared.grid import LatLonGridDescriptor

from ..shared.plot.plotting import plot_polar_comparison, \
    setup_colormap

from ..shared.io.utility import build_config_full_path
from ..shared.io import write_netcdf

from .sea_ice_analysis_task import SeaIceAnalysisTask


class ClimatologyMapSeaIce(SeaIceAnalysisTask):
    """
    General comparison of 2-d model fields against data.  Currently only
    supports sea ice concentration and sea ice thickness

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani
    """

    @classmethod
    def create_tasks(cls, config):  # {{{
        """
        For each comparison grid, construct one task for computing the
        climatologies and one plotting task for each season.  The climatology
        task is a prerequisite of the plotting tasks, but the plotting tasks
        can run in parallel with one another.

        Parameters
        ----------
        config : MpasAnalysisConfigParser object
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """

        # get the name of the class (something like ClimatologyMapSST)
        className = cls.__name__
        sectionPrefix = className[0].lower() + className[1:]

        tasks = []

        # separate climatology tasks for each hemisphere
        for hemisphere in ['NH', 'SH']:
            # convert the first letter to lowercase to get the corresponding
            # seciton name in the config file
            sectionName = '{}{}'.format(sectionPrefix, hemisphere)

            # read in what seasons we want to plot
            seasons = config.getExpression(sectionName, 'seasons')

            if len(seasons) == 0:
                raise ValueError('config section {} does not contain valid '
                                 'list of seasons'.format(sectionName))

            comparisonGridNames = config.getExpression(sectionName,
                                                       'comparisonGrids')

            if len(comparisonGridNames) == 0:
                raise ValueError('config section {} does not contain valid '
                                 'list of comparison grids'.format(
                                         sectionName))

            observationPrefixes = \
                config.getExpression(sectionName, 'observationPrefixes')

            climatologyMapTasks = []
            for comparisonGridName in comparisonGridNames:
                for season in seasons:
                    for observationPrefix in observationPrefixes:
                        # one plotting task for each season, comparison grid,
                        # and type of observations
                        climatologyMapTasks.append(cls(
                                config=config,
                                hemisphere=hemisphere,
                                season=season,
                                comparisonGridName=comparisonGridName,
                                observationPrefix=observationPrefix))

            # we'll use the first task (corresponding to the first season in
            # the list) as a way to get some information we need to build the
            # MpasClimatology task for computing climatologies:
            # mpasVariableName, componentName, fieldName
            firstTask = climatologyMapTasks[0]

            taskSuffix = firstTask.fieldName[0].upper() + \
                firstTask.fieldName[1:] + hemisphere

            climatologyTask = \
                MpasClimatology(config=config,
                                variableList=[firstTask.mpasFieldName],
                                taskSuffix=taskSuffix,
                                componentName=firstTask.componentName,
                                comparisonGridNames=comparisonGridNames,
                                seasons=seasons,
                                tags=['climatology'],
                                iselValues=firstTask.iselValues)

            for index in range(len(climatologyMapTasks)):
                climatologyMapTask = climatologyMapTasks[index]
                # add climatologyTask as a prerequisite of each task so
                # plotting won't happen until we have the required
                # climatologies
                if climatologyMapTask.prerequisiteTasks is None:
                    climatologyMapTask.prerequisiteTasks = \
                        [climatologyTask.taskName]
                else:
                    climatologyMapTask.prerequisiteTasks.append(
                            climatologyTask.taskName)
                # We want to have access to some information from the
                # climatologyTask (namely, we need a way to find out what the
                # names of the climatology files are that it created), so we'll
                # keep a reference to it handy.
                climatologyMapTask.climatologyTask = climatologyTask

            tasks.append(climatologyTask)
            tasks.extend(climatologyMapTasks)

        return tasks  # }}}

    def run(self):  # {{{
        """
        Performs analysis of sea-ice properties by comparing with
        previous model results and/or observations.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of {} climatologies for season {}" \
              "...".format(self.fieldNameInTitle, self.season)

        climatologyTask = self.climatologyTask

        # call setup_and_check() on the climatology task because it will make
        # sure the start and end year are set and correct.  (In parallel mode,
        # this copy of the climatologyTask is different from the one where
        # setup_and_check was already called, and run completed successfully.)
        climatologyTask.setup_and_check()

        self._do_ploting()  # }}}

    def get_task_name_and_tags(self):
        className = type(self).__name__
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        taskPrefix = className[0].lower() + className[1:]

        upperComparison = self.comparisonGridName[0].upper() + \
            self.comparisonGridName[1:]

        taskName = '{}{}{}{}_{}'.format(taskPrefix, self.hemisphere,
                                        self.observationPrefix,
                                        upperComparison, self.season)

        tags = ['climatology', 'horizontalMap', self.comparisonGridName,
                taskPrefix]

        return taskName, tags

    def get_section_name(self):
        className = type(self).__name__
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        sectionName = className[0].lower() + className[1:]

        sectionName = '{}{}'.format(sectionName, self.hemisphere)

        return sectionName

    def _do_ploting(self):  # {{{
        '''
        plots model results, observations and biases.

        Authors
        -------
        Xylar Asay-Davis, Milena Veneziani
        '''
        config = self.config

        mainRunName = config.get('runs', 'mainRunName')

        comparisonDescriptor = \
            get_comparison_descriptor(config, self.comparisonGridName)

        remappedFileName = self.climatologyTask.get_ncclimo_file_name(
                    self.season, 'remapped', comparisonDescriptor)

        obsDescriptor = LatLonGridDescriptor.read(fileName=self.obsFileName,
                                                  latVarName='t_lat',
                                                  lonVarName='t_lon')

        remappedClimatology = xr.open_dataset(remappedFileName)

        modelOutput = remappedClimatology[self.mpasFieldName].values
        if self.maskValue is not None:
            modelOutput = ma.masked_values(modelOutput, self.maskValue)
        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        if self.hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        sectionName = self.get_section_name()

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, sectionName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, sectionName, suffix='Difference')

        referenceLongitude = config.getfloat(sectionName, 'referenceLongitude')
        minimumLatitude = config.getfloat(sectionName, 'minimumLatitude')

        if not os.path.isfile(self.obsFileName):
            raise OSError('Obs file {} not found.'.format(
                self.obsFileName))

        obsRemapper = get_remapper(
                config=config,
                sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(self.fieldName),
                method=config.get('seaIceObservations',
                                  'interpolationMethod'))

        (obsClimatologyFileName, obsRemappedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=self.obsFieldName,
                monthNames=self.season, componentName=self.componentName,
                remapper=obsRemapper)

        if not os.path.exists(obsRemappedFileName):

            # load the observations the first time
            seasonalClimatology = self._build_observational_dataset()
            write_netcdf(seasonalClimatology, obsClimatologyFileName)

            if obsRemapper is None:
                # no need to remap because the observations are on the
                # comparison grid already
                remappedClimatology = seasonalClimatology
            else:
                remappedClimatology = \
                    remap_and_write_climatology(
                        config, seasonalClimatology, obsClimatologyFileName,
                        obsRemappedFileName, obsRemapper)

        else:

            remappedClimatology = xr.open_dataset(obsRemappedFileName)

        observations = remappedClimatology[self.obsFieldName].values
        if self.maskValue is not None:
            observations = ma.masked_values(observations, self.maskValue)

        difference = modelOutput - observations

        startYear = self.climatologyTask.startYear
        endYear = self.climatologyTask.endYear

        title = '{} ({}, years {:04d}-{:04d})'.format(
            self.fieldNameInTitle, self.season, startYear, endYear)
        fileout = '{}/iceconc{}{}_{}_{}_years{:04d}-{:04d}.png'.format(
            self.plotsDirectory, self.observationPrefix, self.hemisphere,
            mainRunName, self.season, startYear, endYear)
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
            obsTitle=self.observationTitleLabel,
            diffTitle='Model-Observations',
            cbarlabel=self.unitsLabel)
        # }}}


class ClimatologyMapSeaIceConc(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice concentration against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, hemisphere, season, comparisonGridName,
                 observationPrefix):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        seasons : ['NH', 'SH']
            The hemisphere to plot

        seasons : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        observationPrefix : {'NASATeam', 'Bootstrap'}
            A prefix describing the set of observations to use

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'seaIceConc'
        self.fieldNameInTitle = 'Sea ice concentration'
        self.hemisphere = hemisphere
        self.season = season
        self.comparisonGridName = comparisonGridName
        self.observationPrefix = observationPrefix

        self.mpasFieldName = 'timeMonthly_avg_iceAreaCell'
        self.iselValues = None

        taskName, tags = self.get_task_name_and_tags()

        tags.append(self.fieldName)

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceConc, self).__init__(config=config,
                                                       taskName=taskName,
                                                       componentName='seaIce',
                                                       tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceConc, self).setup_and_check()

        self.obsFileName = build_config_full_path(
            self.config, 'seaIceObservations',
            'concentration{}{}_{}'.format(self.observationPrefix,
                                          self.hemisphere,
                                          self.season))

        self.obsFieldName = 'AICE'

        self.observationTitleLabel = \
            'Observations (SSM/I {})'.format(self.observationPrefix)

        self.outFileLabel = 'iceconc{}{}'.format(self.observationPrefix,
                                                 self.hemisphere)
        self.unitsLabel = 'fraction'

        self.maskValue = None

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(self.obsFileName)

        return dsObs  # }}}

    # }}}


class ClimatologyMapSeaIceThick(ClimatologyMapSeaIce):  # {{{
    """
    An analysis task for comparison of sea ice thickness against observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, hemisphere, season, comparisonGridName,
                 observationPrefix=''):
        # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        seasons : ['NH', 'SH']
            The hemisphere to plot

        seasons : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        observationPrefix : {''}
            A prefix describing the set of observations to use (empty for
            this task)

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'seaIceThick'
        self.fieldNameInTitle = 'Sea ice thickness'
        self.hemisphere = hemisphere
        self.season = season
        self.comparisonGridName = comparisonGridName
        self.observationPrefix = observationPrefix

        self.mpasFieldName = 'timeMonthly_avg_iceVolumeCell'
        self.iselValues = None

        taskName, tags = self.get_task_name_and_tags()

        tags.append(self.fieldName)

        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSeaIceThick, self).__init__(config=config,
                                                        taskName=taskName,
                                                        componentName='seaIce',
                                                        tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class
        # (ClimatologyMapSeaIce), which will perform some common setup,
        # including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSeaIceThick, self).setup_and_check()

        self.obsFileName = build_config_full_path(
            self.config, 'seaIceObservations',
            'thickness{}_{}'.format(self.hemisphere, self.season))

        self.obsFieldName = 'HI'

        self.observationTitleLabel = 'Observations (ICESat)'

        self.outFileLabel = 'icethick{}'.format(self.hemisphere)

        self.unitsLabel = 'm'

        self.maskValue = 0

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        dsObs = xr.open_mfdataset(self.obsFileName)

        return dsObs  # }}}

    # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
