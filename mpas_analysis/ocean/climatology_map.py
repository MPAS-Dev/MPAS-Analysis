"""
An analysis tasks for comparison of 2D model fields against observations.
Currently supports sea surface temperature (sst), sea surface salinity
(sss) and mixed layer depth (mld)

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
"""

from __future__ import absolute_import, division, print_function, unicode_literals

import xarray as xr
import datetime
import numpy as np
import os

from ..shared import AnalysisTask

from ..shared.plot.plotting import plot_global_comparison, \
    setup_colormap
from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path, make_directories
from ..shared.html import write_image_xml

from ..shared.climatology import get_comparison_descriptor, \
    get_remapper, get_observation_climatology_file_names, \
    compute_climatology, remap_and_write_climatology, \
    RemapMpasClimatologySubtask

from ..shared.grid import LatLonGridDescriptor

from ..shared.mpas_xarray import mpas_xarray

from ..shared.interpolation import Remapper


class ClimatologyMapOcean(AnalysisTask):  # {{{
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

    galleryGroup : str
        In the generated website, the name of the group of analysis galleries
        to which this analysis belongs

    groupLink : str
        the link within the component webpage used to identify the gallery
        group (a short version of the gallery group name with no spaces or
        punctuation)

    groupSubtitle : str
        a subtitle for the gallery group

    galleryName : str, optional
        the name of the gallery (or possibly a subtitle for the gallery group
        if there is only one gallery in the group)

    remapClimatologySubtask : ``RemapMpasClimatologySubtask``
        The subtask that remaps the climatologies this task will plot.
        The ``remapClimatologySubtask`` is needed to determine the file names
        of the climatology output.

    mpasClimatologyTask : ``MpasClimatologyTask``
        The task that produced the climatology to be remapped and plotted

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
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
        super(ClimatologyMapOcean, self).__init__(config=config,
                                                  taskName=taskName,
                                                  componentName='ocean',
                                                  tags=tags)

        self.mpasClimatologyTask = mpasClimatologyTask

        sectionName = self.taskName

        # read in what seasons we want to plot
        seasons = config.getExpression(sectionName, 'seasons')

        if len(seasons) == 0:
            raise ValueError('config section {} does not contain valid list '
                             'of seasons'.format(sectionName))

        # comparisonGridNames = config.getExpression(sectionName,
        #                                           'comparisonGrids')

        # if len(comparisonGridNames) == 0:
        #     raise ValueError('config section {} does not contain valid list '
        #                      'of comparison grids'.format(sectionName))

        # the variable self.mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        self.remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=self.mpasClimatologyTask,
            parentTask=self,
            climatologyName=self.fieldName,
            variableList=[self.mpasFieldName],
            # comparisonGridNames=comparisonGridNames,
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
        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapOcean, self).setup_and_check()

        config = self.config
        self.startYear = self.mpasClimatologyTask.startYear
        self.startDate = self.mpasClimatologyTask.startDate
        self.endYear = self.mpasClimatologyTask.endYear
        self.endDate = self.mpasClimatologyTask.endDate

        mainRunName = config.get('runs', 'mainRunName')
        seasons = config.getExpression(self.taskName, 'seasons')

        # we set up the remapper here because ESFM_RegridWeightGen seems to
        # have trouble if it runs in another process (or in several at once)
        self._setup_obs_remapper()

        self.xmlFileNames = []
        self.filePrefixes = {}

        for season in seasons:
            filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                    self.outFileLabel, mainRunName,
                    season, self.startYear, self.endYear)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[season] = filePrefix

        # make the mapping directory, because doing so within each process
        # seems to be giving ESMF_RegridWeightGen some trouble
        mappingSubdirectory = build_config_full_path(config, 'output',
                                                     'mappingSubdirectory')
        make_directories(mappingSubdirectory)

        # }}}

    def run_task(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
        """

        self.logger.info("\nPlotting 2-d maps of {} climatologies...".format(
            self.fieldNameInTitle))

        # get local versions of member variables for convenience
        config = self.config
        fieldName = self.fieldName

        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for ocn_modelvsobs calculation')

        mainRunName = self.config.get('runs', 'mainRunName')

        seasons = config.getExpression(self.taskName, 'seasons')

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, self.taskName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, self.taskName, suffix='Difference')

        dsRestart = xr.open_dataset(restartFileName)
        dsRestart = mpas_xarray.subset_variables(dsRestart, ['maxLevelCell'])

        dsObs = None

        # Interpolate and compute biases
        for season in seasons:

            monthValues = constants.monthDictionary[season]

            remappedFileName = self.remapClimatologySubtask.get_file_name(
                season=season, stage='remapped', comparisonGridName='latlon')

            remappedClimatology = xr.open_dataset(remappedFileName)

            modelOutput = \
                remappedClimatology[self.mpasFieldName].values

            lon = remappedClimatology['lon'].values
            lat = remappedClimatology['lat'].values

            lonTarg, latTarg = np.meshgrid(lon, lat)

            # now the observations
            (climatologyFileName, remappedFileName) = \
                get_observation_climatology_file_names(
                    config=config, fieldName=fieldName, monthNames=season,
                    componentName='ocean', remapper=self.obsRemapper)

            if not os.path.exists(remappedFileName):

                if dsObs is None:
                    # load the observations the first time
                    dsObs = self._build_observational_dataset()

                seasonalClimatology = compute_climatology(
                    dsObs, monthValues, maskVaries=True)

                if self.obsRemapper is None:
                    # no need to remap because the observations are on the
                    # comparison grid already
                    remappedClimatology = seasonalClimatology
                else:
                    remappedClimatology = \
                        remap_and_write_climatology(
                            config, seasonalClimatology, climatologyFileName,
                            remappedFileName, self.obsRemapper,
                            logger=self.logger)

            else:

                remappedClimatology = xr.open_dataset(remappedFileName)
            observations = remappedClimatology[self.obsFieldName].values

            bias = modelOutput - observations

            filePrefix = self.filePrefixes[season]
            outFileName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)
            title = '{} ({}, years {:04d}-{:04d})'.format(
                    self.fieldNameInTitle, season, self.startYear,
                    self.endYear)
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

            caption = '{} {}'.format(season, self.imageCaption)
            write_image_xml(
                config,
                filePrefix,
                componentName='Ocean',
                componentSubdirectory='ocean',
                galleryGroup=self.galleryGroup,
                groupSubtitle=self.groupSubtitle,
                groupLink=self.groupLink,
                gallery=self.galleryName,
                thumbnailDescription=season,
                imageDescription=caption,
                imageCaption=caption)

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
        fieldName = self.fieldName

        seasons = config.getExpression(self.taskName, 'seasons')

        # make reamppers
        comparisonDescriptor = get_comparison_descriptor(
                config=config, comparisonGridName='latlon')
        self.comparisonGridName = comparisonDescriptor.meshName

        obsDescriptor = LatLonGridDescriptor.read(fileName=self.obsFileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')

        origObsRemapper = Remapper(comparisonDescriptor, obsDescriptor)

        season = seasons[0]

        # now the observations
        (climatologyFileName, remappedFileName) = \
            get_observation_climatology_file_names(
                config=config, fieldName=fieldName, monthNames=season,
                componentName='ocean', remapper=origObsRemapper)

        # make the remapper for the climatology
        self.obsRemapper = get_remapper(
            config=config, sourceDescriptor=obsDescriptor,
            comparisonDescriptor=comparisonDescriptor,
            mappingFilePrefix='map_obs_{}'.format(fieldName),
            method=config.get('oceanObservations',
                              'interpolationMethod'),
            logger=self.logger)

        # }}}
    # }}}


class ClimatologyMapSST(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface temperature (sst) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
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
        self.fieldName = 'sst'
        self.fieldNameInTitle = 'SST'
        self.mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        self.iselValues = {'nVertLevels': 0}

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSST, self).__init__(
            config=config,
            mpasClimatologyTask=mpasClimatologyTask,
            taskName='climatologyMapSST',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """

        self.outFileLabel = 'sstHADOI'

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
        self.unitsLabel = r'$^o$C'

        # variables for XML and webpages
        self.galleryGroup = 'Global Sea Surface Temperature'
        self.groupSubtitle = None
        self.groupLink = 'sst'
        self.galleryName = 'Observations: Hadley-NOAA-OI'
        self.imageCaption = 'Mean Sea Surface Temperature'

        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapSST, self).setup_and_check()

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
        self.fieldName = 'sss'
        self.fieldNameInTitle = 'SSS'
        self.mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        self.iselValues = {'nVertLevels': 0}

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSSS, self).__init__(
            config=config,
            mpasClimatologyTask=mpasClimatologyTask,
            taskName='climatologyMapSSS',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """

        self.outFileLabel = 'sssAquarius'

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/Aquarius_V3_SSS_Monthly.nc'.format(
                observationsDirectory)

        self.obsFieldName = 'SSS'

        self.observationTitleLabel = 'Observations (Aquarius, 2011-2014)'
        self.unitsLabel = 'PSU'

        # variables for XML and webpages
        self.galleryGroup = 'Global Sea Surface Salinity'
        self.groupSubtitle = None
        self.groupLink = 'sss'
        self.galleryName = 'Observations: Aquarius'
        self.imageCaption = 'Mean Sea Surface Salinity'

        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapSSS, self).setup_and_check()

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

        self.fieldName = 'mld'
        self.fieldNameInTitle = 'MLD'
        self.mpasFieldName = 'timeMonthly_avg_dThreshMLD'
        self.iselValues = None

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapMLD, self).__init__(
            config=config,
            mpasClimatologyTask=mpasClimatologyTask,
            taskName='climatologyMapMLD',
            tags=['climatology', 'horizontalMap', self.fieldName])

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        """

        self.outFileLabel = 'mldHolteTalleyARGO'

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
        self.unitsLabel = 'm'

        # variables for XML and webpages
        self.galleryGroup = 'Global Mixed-Layer Depth'
        self.groupSubtitle = None
        self.groupLink = 'mld'
        self.galleryName = 'Observations: Holte-Talley ARGO'
        self.imageCaption = 'Mean Mixed-Layer Depth'

        # call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapMLD, self).setup_and_check()

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
