"""
An analysis tasks for comparison of 2D model fields against observations.
Currently supports sea surface temperature (sst), sea surface salinity
(sss) and mixed layer depth (mld)

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
"""

import xarray as xr
import datetime
import numpy as np
import os

from ..shared.analysis_task import AnalysisTask

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap
from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path
from ..shared.io import write_netcdf

from ..shared.climatology import get_comparison_descriptor, \
    get_remapper, get_observation_climatology_file_names, \
    compute_climatology, remap_and_write_climatology, MpasClimatology

from ..shared.grid import LatLonGridDescriptor

from ..shared.mpas_xarray import mpas_xarray

from ..shared.interpolation import Remapper


class ClimatologyMapOcean(AnalysisTask):  # {{{
    """
    An analysis task for comparison of 2D model fields against observations.

    Attributes
    ----------
    climatologyTask : MpasClimatology object
        The climatology task that has already finished by the time the ``run``
        method gets called in this class.  It is used to figure out the
        names of the climatology files that it created.

    seasons : str
        A season (keys in ``shared.constants.monthDictionary``) to be plotted

    comparisonGridName : {'latlon', 'antarctic'}
        The name of the comparison grid that where model data an observations
        are compared

    fieldName : str
        A short name of the field being analyzed

    fieldNameInTitle : str
        An equivalent name of the field appropriate for figure titles

    mpasFieldName : str
        The name of the MPAS timeSeriesStatsMonthly variable to be analyzed

    iselValues : dict
        A dictionary of dimensions and indices (or ``None``) used to extract
        a slice of the MPAS field.

    obsFileName : str
        A file containing observtions from which to construct seasonal
        climatologies.

    obsFieldName : str
        A variable in the observations file to plot

    observationTitleLabel : str
        A label on the subplot showing observations

    outFileLabel : str
        A prefix for the resulting image file.

    unitsLabel : str
        The units for the field being plotted.

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
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
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        sectionName = className[0].lower() + className[1:]

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

        climatologyMapTasks = []
        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # one plotting task for each season and comparison grid
                climatologyMapTasks.append(cls(
                        config=config, season=season,
                        comparisonGridName=comparisonGridName))

        # we'll use the first task (corresponding to the first season in
        # the list) as a way to get some information we need to build the
        # MpasClimatology task for computing climatologies:
        # mpasVariableName, componentName, fieldName
        firstTask = climatologyMapTasks[0]

        taskSuffix = firstTask.fieldName[0].upper() + \
            firstTask.fieldName[1:]

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

        tasks = [climatologyTask] + climatologyMapTasks

        return tasks  # }}}

    def run(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of {} climatologies for season {}" \
              "...".format(self.fieldNameInTitle, self.season)

        # get local versions of member variables for convenience
        config = self.config
        fieldName = self.fieldName
        climatologyTask = self.climatologyTask

        # call setup_and_check() on the climatology task because it will make
        # sure the start and end year are set and correct.  (In parallel mode,
        # this copy of the climatologyTask is different from the one where
        # setup_and_check was already called, and run completed successfully.)
        climatologyTask.setup_and_check()

        mainRunName = config.get('runs', 'mainRunName')

        obsDescriptor = LatLonGridDescriptor.read(fileName=self.obsFileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')

        comparisonDescriptor = \
            get_comparison_descriptor(config, self.comparisonGridName)

        origObsRemapper = Remapper(comparisonDescriptor, obsDescriptor)

        sectionName = self.get_section_name()

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, sectionName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, sectionName, suffix='Difference')

        startYear = climatologyTask.startYear
        endYear = climatologyTask.endYear

        # Interpolate and compute biases
        monthValues = constants.monthDictionary[self.season]

        remappedFileName = climatologyTask.get_ncclimo_file_name(
                    self.season, 'remapped', comparisonDescriptor)

        remappedClimatology = xr.open_dataset(remappedFileName)

        modelOutput = \
            remappedClimatology[self.mpasFieldName].values

        lon = remappedClimatology['lon'].values
        lat = remappedClimatology['lat'].values

        lonTarg, latTarg = np.meshgrid(lon, lat)

        # now the observations
        (climatologyFileName, remappedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=fieldName, monthNames=self.season,
                componentName='ocean', remapper=origObsRemapper)

        if not os.path.exists(remappedFileName):

            # load the observations the first time
            dsObs = self._build_observational_dataset()

            seasonalClimatology = compute_climatology(
                dsObs, monthValues, maskVaries=True)

            seasonalClimatology.load()
            seasonalClimatology.close()
            write_netcdf(seasonalClimatology, climatologyFileName)
            # make the remapper for the climatology
            obsDescriptor = LatLonGridDescriptor.read(
                    fileName=climatologyFileName,
                    latVarName='lat',
                    lonVarName='lon')

            obsRemapper = get_remapper(
                config=config, sourceDescriptor=obsDescriptor,
                comparisonDescriptor=comparisonDescriptor,
                mappingFilePrefix='map_obs_{}'.format(fieldName),
                method=config.get('oceanObservations',
                                  'interpolationMethod'))

            if obsRemapper is None:
                # no need to remap because the observations are on the
                # comparison grid already
                remappedClimatology = seasonalClimatology
            else:
                remappedClimatology = \
                    remap_and_write_climatology(
                        config, seasonalClimatology, climatologyFileName,
                        remappedFileName, obsRemapper)

        else:

            remappedClimatology = xr.open_dataset(remappedFileName)
        observations = remappedClimatology[self.obsFieldName].values

        bias = modelOutput - observations

        outFileName = '{}/{}_{}_{}_years{:04d}-{:04d}.png'.format(
                self.plotsDirectory, self.outFileLabel, mainRunName,
                self.season, startYear, endYear)
        title = '{} ({}, years {:04d}-{:04d})'.format(
                self.fieldNameInTitle, self.season, startYear, endYear)
        plot_global_comparison(config,
                               lonTarg,
                               latTarg,
                               modelOutput,
                               observations,
                               bias,
                               colormapResult,
                               colorbarLevelsResult,
                               colormapDifference,
                               colorbarLevelsDifference,
                               fileout=outFileName,
                               title=title,
                               modelTitle='{}'.format(mainRunName),
                               obsTitle=self.observationTitleLabel,
                               diffTitle='Model-Observations',
                               cbarlabel=self.unitsLabel)

        # }}}

    def get_task_name_and_tags(self):
        className = type(self).__name__
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        taskPrefix = className[0].lower() + className[1:]

        upperComparison = self.comparisonGridName[0].upper() + \
            self.comparisonGridName[1:]

        taskName = '{}{}_{}'.format(taskPrefix, upperComparison, self.season)

        tags = ['climatology', 'horizontalMap', self.comparisonGridName,
                taskPrefix]

        return taskName, tags

    def get_section_name(self):
        className = type(self).__name__
        # convert the first letter to lowercase to get the corresponding
        # seciton name in the config file
        sectionName = className[0].lower() + className[1:]

        upperComparison = self.comparisonGridName[0].upper() + \
            self.comparisonGridName[1:]

        comparisonSectionName = '{}{}'.format(sectionName, upperComparison)

        if self.config.has_section(comparisonSectionName):
            sectionName = comparisonSectionName

        return sectionName

    # }}}


class ClimatologyMapSST(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface temperature (sst) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, season, comparisonGridName):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        seasons : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sst'
        self.fieldNameInTitle = 'SST'
        self.season = season
        self.comparisonGridName = comparisonGridName

        self.mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        self.iselValues = {'nVertLevels': 0}

        taskName, tags = self.get_task_name_and_tags()

        tags.append(self.fieldName)

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSST, self).__init__(config=config,
                                                taskName=taskName,
                                                componentName='ocean',
                                                tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSST, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(
                observationsDirectory)

        self.obsFieldName = 'SST'

        climStartYear = self.config.getint('oceanObservations',
                                           'sstClimatologyStartYear')
        climEndYear = self.config.getint('oceanObservations',
                                         'sstClimatologyEndYear')

        if climStartYear < 1925:
            period = 'pre-industrial'
        else:
            period = 'present-day'

        # Set appropriate figure labels for SST
        self.observationTitleLabel = \
            'Observations (Hadley/OI, {} {:04d}-{:04d})'.format(period,
                                                                climStartYear,
                                                                climEndYear)
        self.outFileLabel = 'sstHADOI'
        self.unitsLabel = r'$^o$C'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        climStartYear = self.config.getint('oceanObservations',
                                           'sstClimatologyStartYear')
        climEndYear = self.config.getint('oceanObservations',
                                         'sstClimatologyEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_mfdataset(self.obsFileName)
        dsObs.rename({'time': 'Time'}, inplace=True)
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs  # }}}

    # }}}


class ClimatologyMapSSS(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface salinity (sss) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, season, comparisonGridName):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        seasons : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sss'
        self.fieldNameInTitle = 'SSS'
        self.season = season
        self.comparisonGridName = comparisonGridName

        self.mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        self.iselValues = {'nVertLevels': 0}

        taskName, tags = self.get_task_name_and_tags()

        tags.append(self.fieldName)

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSSS, self).__init__(config=config,
                                                taskName=taskName,
                                                componentName='ocean',
                                                tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapSSS, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/Aquarius_V3_SSS_Monthly.nc'.format(
                observationsDirectory)

        self.obsFieldName = 'SSS'

        self.observationTitleLabel = 'Observations (Aquarius, 2011-2014)'
        self.outFileLabel = 'sssAquarius'
        self.unitsLabel = 'PSU'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        timeStart = datetime.datetime(2011, 8, 1)
        timeEnd = datetime.datetime(2014, 12, 31)

        dsObs = xr.open_mfdataset(self.obsFileName)
        dsObs.rename({'time': 'Time'}, inplace=True)
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs  # }}}

    # }}}


class ClimatologyMapMLD(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of mixed layer depth (mld) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config, season, comparisonGridName):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        seasons : str
            A season (keys in ``shared.constants.monthDictionary``) to be
            plotted

        comparisonGridName : {'latlon', 'antarctic'}
            The name of the comparison grid that where model data an
            observations are compared

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'mld'
        self.fieldNameInTitle = 'MLD'
        self.season = season
        self.comparisonGridName = comparisonGridName

        self.mpasFieldName = 'timeMonthly_avg_dThreshMLD'
        self.iselValues = None

        taskName, tags = self.get_task_name_and_tags()

        tags.append(self.fieldName)

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapMLD, self).__init__(config=config,
                                                taskName=taskName,
                                                componentName='ocean',
                                                tags=tags)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(ClimatologyMapMLD, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/holtetalley_mld_climatology.nc'.format(
                observationsDirectory)

        self.obsFieldName = 'mld_dt_mean'

        # Set appropriate MLD figure labels
        self.observationTitleLabel = \
            'Observations (HolteTalley density threshold MLD)'
        self.outFileLabel = 'mldHolteTalleyARGO'
        self.unitsLabel = 'm'

        # }}}

    def _build_observational_dataset(self):  # {{{
        '''
        read in the data sets for observations, and possibly rename some
        variables and dimensions

        Authors
        -------
        Xylar Asay-Davis
        '''

        # Load MLD observational data
        dsObs = xr.open_mfdataset(self.obsFileName)

        # Increment month value to be consistent with the model output
        dsObs.iMONTH.values += 1
        # Rename the dimensions to be consistent with other obs. data sets
        dsObs.rename({'month': 'calmonth', 'lat': 'latCoord',
                      'lon': 'lonCoord'}, inplace=True)
        dsObs.rename({'iMONTH': 'Time', 'iLAT': 'lat', 'iLON': 'lon'},
                     inplace=True)

        # set the coordinates now that the dimensions have the same names
        dsObs.coords['lat'] = dsObs['latCoord']
        dsObs.coords['lon'] = dsObs['lonCoord']
        dsObs.coords['Time'] = dsObs['calmonth']
        dsObs.coords['month'] = ('Time', np.array(dsObs['calmonth'], int))

        # no meaningful year since this is already a climatology
        dsObs.coords['year'] = ('Time', np.ones(dsObs.dims['Time'], int))

        dsObs = mpas_xarray.subset_variables(dsObs, [self.obsFieldName,
                                                     'month'])

        return dsObs  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
