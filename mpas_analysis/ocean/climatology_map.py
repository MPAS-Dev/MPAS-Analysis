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
from ..shared.html import write_image_xml

from ..shared.climatology import get_lat_lon_comparison_descriptor, \
    get_remapper, get_mpas_climatology_dir_name, \
    get_observation_climatology_file_names, compute_climatology, \
    remap_and_write_climatology, update_climatology_bounds_from_file_names, \
    compute_climatologies_with_ncclimo, get_ncclimo_season_file_name

from ..shared.grid import MpasMeshDescriptor, LatLonGridDescriptor

from ..shared.mpas_xarray import mpas_xarray

from ..shared.interpolation import Remapper


class ClimatologyMapOcean(AnalysisTask):  # {{{
    """
    An analysis task for comparison of 2D model fields against observations.

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """

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
        #     self.calendar
        super(ClimatologyMapOcean, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = self.config.get('climatology', 'startDate')
        endDate = self.config.get('climatology', 'endDate')
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.inputFiles = self.historyStreams.readpath(
                streamName, startDate=startDate, endDate=endDate,
                calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, startDate, endDate))

        changed, self.startYear, self.endYear, self.startDate, self.endDate = \
            update_climatology_bounds_from_file_names(self.inputFiles,
                                                      self.config)
        mainRunName = self.config.get('runs', 'mainRunName')
        comparisonTimes = self.config.getExpression(self.taskName,
                                                    'comparisonTimes')

        self.xmlFileNames = []
        self.filePrefixes = {}

        for season in comparisonTimes:
            filePrefix = '{}_{}_{}_years{:04d}-{:04d}'.format(
                    self.outFileLabel, mainRunName,
                    season, self.startYear, self.endYear)
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))
            self.filePrefixes[season] = filePrefix

        # }}}

    def run(self):  # {{{
        """
        Plots a comparison of ACME/MPAS output to SST, MLD or SSS observations

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
        """

        print "\nPlotting 2-d maps of {} climatologies...".format(
            self.fieldNameInTitle)

        # get local versions of member variables for convenience
        config = self.config
        fieldName = self.fieldName

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

        mainRunName = config.get('runs', 'mainRunName')

        try:
            restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for ocn_modelvsobs calculation')

        outputTimes = config.getExpression(self.taskName, 'comparisonTimes')

        comparisonDescriptor = get_lat_lon_comparison_descriptor(config)

        mpasDescriptor = MpasMeshDescriptor(
            restartFileName, meshName=config.get('input', 'mpasMeshName'))

        parallel = self.config.getint('execute', 'parallelTaskCount') > 1
        if parallel:
            # avoid writing the same mapping file from multiple processes
            mappingFilePrefix = 'map_{}'.format(self.taskName)
        else:
            mappingFilePrefix = 'map'

        mpasRemapper = get_remapper(
            config=config, sourceDescriptor=mpasDescriptor,
            comparisonDescriptor=comparisonDescriptor,
            mappingFilePrefix=mappingFilePrefix,
            method=config.get('climatology', 'mpasInterpolationMethod'))

        obsDescriptor = LatLonGridDescriptor.read(fileName=self.obsFileName,
                                                  latVarName='lat',
                                                  lonVarName='lon')

        origObsRemapper = Remapper(comparisonDescriptor, obsDescriptor)

        (colormapResult, colorbarLevelsResult) = setup_colormap(
            config, self.taskName, suffix='Result')
        (colormapDifference, colorbarLevelsDifference) = setup_colormap(
            config, self.taskName, suffix='Difference')

        climatologyDirectory = \
            get_mpas_climatology_dir_name(
                config=config,
                fieldName=self.mpasFieldName,
                mpasMeshName=mpasDescriptor.meshName)

        (maskedClimatologyDirectory, remappedDirectory) = \
            get_mpas_climatology_dir_name(
                config=config,
                fieldName='{}_masked'.format(fieldName),
                mpasMeshName=mpasDescriptor.meshName,
                comparisonGridName=comparisonDescriptor.meshName)

        dsRestart = xr.open_dataset(restartFileName)
        dsRestart = mpas_xarray.subset_variables(dsRestart, ['maxLevelCell'])

        startYear = config.getint('climatology', 'startYear')
        endYear = config.getint('climatology', 'endYear')

        modelName = 'mpaso'

        allExist = True
        for season in outputTimes:
            climatologyFileName = get_ncclimo_season_file_name(
                    climatologyDirectory, modelName, season, self.startYear,
                    self.endYear)
            if not os.path.exists(climatologyFileName):
                allExist = False
                break

        if not allExist:

            compute_climatologies_with_ncclimo(
                    config=config,
                    inDirectory=self.historyDirectory,
                    outDirectory=climatologyDirectory,
                    startYear=startYear,
                    endYear=endYear,
                    variableList=[self.mpasFieldName],
                    modelName=modelName,
                    seasons=outputTimes,
                    decemberMode='sdd')

        dsObs = None
        obsRemapperBuilt = False

        # Interpolate and compute biases
        for season in outputTimes:

            monthValues = constants.monthDictionary[season]

            climatologyFileName = \
                get_ncclimo_season_file_name(climatologyDirectory, modelName,
                                             season, startYear, endYear)

            maskedClimatologyFileName = \
                get_ncclimo_season_file_name(maskedClimatologyDirectory,
                                             modelName, season, startYear,
                                             endYear)

            remappedFileName = \
                get_ncclimo_season_file_name(remappedDirectory,
                                             modelName, season, startYear,
                                             endYear)

            if not os.path.exists(maskedClimatologyFileName):
                # slice and mask the data set
                climatology = xr.open_dataset(climatologyFileName)
                iselValues = {'Time': 0}
                if self.iselValues is not None:
                    iselValues.update(self.iselValues)
                # select only Time=0 and possibly only the desired vertical
                # slice
                climatology = climatology.isel(**iselValues)

                # mask the data set
                climatology[self.mpasFieldName] = \
                    climatology[self.mpasFieldName].where(
                            dsRestart.maxLevelCell > 0)

                write_netcdf(climatology, maskedClimatologyFileName)

            if not os.path.exists(remappedFileName):
                mpasRemapper.remap_file(inFileName=maskedClimatologyFileName,
                                        outFileName=remappedFileName,
                                        overwrite=True)

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
                    componentName='ocean', remapper=origObsRemapper)

            if not os.path.exists(remappedFileName):

                if dsObs is None:
                    # load the observations the first time
                    dsObs = self._build_observational_dataset()

                seasonalClimatology = compute_climatology(
                    dsObs, monthValues, maskVaries=True)

                if not obsRemapperBuilt:
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

                    obsRemapperBuilt = True

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

    # }}}


class ClimatologyMapSST(ClimatologyMapOcean):  # {{{
    """
    An analysis task for comparison of sea surface temperature (sst) against
    observations

    Authors
    -------
    Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani
    """
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sst'
        self.fieldNameInTitle = 'SST'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSST, self).__init__(
            config=config,
            taskName='climatologyMapSST',
            componentName='ocean',
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

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapSST, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            "{}/MODEL.SST.HAD187001-198110.OI198111-201203.nc".format(
                observationsDirectory)

        self.mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        self.iselValues = {'nVertLevels': 0}

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
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """
        self.fieldName = 'sss'
        self.fieldNameInTitle = 'SSS'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapSSS, self).__init__(
            config=config,
            taskName='climatologyMapSSS',
            componentName='ocean',
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

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapSSS, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/Aquarius_V3_SSS_Monthly.nc'.format(
                observationsDirectory)

        self.mpasFieldName = 'timeMonthly_avg_activeTracers_salinity'
        self.iselValues = {'nVertLevels': 0}

        self.obsFieldName = 'SSS'

        self.observationTitleLabel = 'Observations (Aquarius, 2011-2014)'
        self.unitsLabel = 'PSU'

        # variables for XML and webpages
        self.galleryGroup = 'Global Sea Surface Salinity'
        self.groupSubtitle = None
        self.groupLink = 'sss'
        self.galleryName = 'Observations: Aquarius'
        self.imageCaption = 'Mean Sea Surface Salinity'

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
    def __init__(self, config):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis
        """

        self.fieldName = 'mld'
        self.fieldNameInTitle = 'MLD'

        # call the constructor from the base class (ClimatologyMapOcean)
        super(ClimatologyMapMLD, self).__init__(
            config=config,
            taskName='climatologyMapMLD',
            componentName='ocean',
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

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(ClimatologyMapMLD, self).setup_and_check()

        observationsDirectory = build_config_full_path(
            self.config, 'oceanObservations',
            '{}Subdirectory'.format(self.fieldName))

        self.obsFileName = \
            '{}/holtetalley_mld_climatology.nc'.format(
                observationsDirectory)

        self.mpasFieldName = 'timeMonthly_avg_dThreshMLD'
        self.iselValues = None

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
