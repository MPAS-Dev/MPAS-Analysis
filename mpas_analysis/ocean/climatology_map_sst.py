# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
import xarray as xr
import datetime

from pyremap import LatLonGridDescriptor

from mpas_analysis.shared import AnalysisTask

from mpas_analysis.shared.io.utility import build_obs_path

from mpas_analysis.shared.climatology import RemapMpasClimatologySubtask, \
    RemapObservedClimatologySubtask

from mpas_analysis.ocean.plot_climatology_map_subtask import \
    PlotClimatologyMapSubtask


class ClimatologyMapSST(AnalysisTask):
    """
    An analysis task for comparison of sea surface temperature (sst) against
    observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

    def __init__(self, config, mpasClimatologyTask,
                 controlConfig=None):
        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasClimatologyTask : ``MpasClimatologyTask``
            The task that produced the climatology to be remapped and plotted

        controlconfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        fieldName = 'sst'
        # call the constructor from the base class (AnalysisTask)
        super(ClimatologyMapSST, self).__init__(
            config=config, taskName='climatologyMapSST',
            componentName='ocean',
            tags=['climatology', 'horizontalMap', fieldName, 'publicObs'])

        mpasFieldName = 'timeMonthly_avg_activeTracers_temperature'
        iselValues = {'nVertLevels': 0}

        sectionName = self.taskName

        climStartYear = config.getint(sectionName, 'obsStartYear')
        climEndYear = config.getint(sectionName, 'obsEndYear')

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

        # the variable mpasFieldName will be added to mpasClimatologyTask
        # along with the seasons.
        remapClimatologySubtask = RemapMpasClimatologySubtask(
            mpasClimatologyTask=mpasClimatologyTask,
            parentTask=self,
            climatologyName=fieldName,
            variableList=[mpasFieldName],
            comparisonGridNames=comparisonGridNames,
            seasons=seasons,
            iselValues=iselValues)

        if controlConfig is None:
            if climStartYear < 1925:
                period = 'pre-industrial'
            else:
                period = 'present-day'

            refTitleLabel = \
                'Observations (Hadley/OI, {} {:04d}-{:04d})'.format(
                    period, climStartYear, climEndYear)

            observationsDirectory = build_obs_path(
                config, 'ocean', '{}Subdirectory'.format(fieldName))

            obsFileName = \
                "{}/MODEL.SST.HAD187001-198110.OI198111-201203_" \
                "20180710.nc".format(observationsDirectory)
            refFieldName = 'sst'
            outFileLabel = 'sstHADOI'
            galleryName = 'Observations: Hadley-NOAA-OI'

            remapObservationsSubtask = RemapObservedSSTClimatology(
                parentTask=self, seasons=seasons, fileName=obsFileName,
                outFilePrefix=refFieldName,
                comparisonGridNames=comparisonGridNames)
            self.add_subtask(remapObservationsSubtask)
            diffTitleLabel = 'Model - Observations'

        else:
            remapObservationsSubtask = None
            controlRunName = controlConfig.get('runs', 'mainRunName')
            galleryName = None
            refTitleLabel = 'Control: {}'.format(controlRunName)

            refFieldName = mpasFieldName
            outFileLabel = 'sst'
            diffTitleLabel = 'Main - Control'

        for comparisonGridName in comparisonGridNames:
            for season in seasons:
                # make a new subtask for this season and comparison grid
                subtask = PlotClimatologyMapSubtask(
                    self, season, comparisonGridName, remapClimatologySubtask,
                    remapObservationsSubtask, controlConfig=controlConfig)

                subtask.set_plot_info(
                    outFileLabel=outFileLabel,
                    fieldNameInTitle='SST',
                    mpasFieldName=mpasFieldName,
                    refFieldName=refFieldName,
                    refTitleLabel=refTitleLabel,
                    diffTitleLabel=diffTitleLabel,
                    unitsLabel=r'$^o$C',
                    imageCaption='Mean Sea Surface Temperature',
                    galleryGroup='Sea Surface Temperature',
                    groupSubtitle=None,
                    groupLink='sst',
                    galleryName=galleryName)

                self.add_subtask(subtask)


class RemapObservedSSTClimatology(RemapObservedClimatologySubtask):
    """
    A subtask for reading and remapping SST observations
    """
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

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
                                                  latVarName='lat',
                                                  lonVarName='lon')
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

        sectionName = self.taskName
        climStartYear = self.config.getint(sectionName, 'obsStartYear')
        climEndYear = self.config.getint(sectionName, 'obsEndYear')
        timeStart = datetime.datetime(year=climStartYear, month=1, day=1)
        timeEnd = datetime.datetime(year=climEndYear, month=12, day=31)

        dsObs = xr.open_dataset(fileName)
        dsObs = dsObs.rename({'time': 'Time', 'SST': 'sst'})
        dsObs = dsObs.sel(Time=slice(timeStart, timeEnd))
        dsObs.coords['month'] = dsObs['Time.month']
        dsObs.coords['year'] = dsObs['Time.year']

        return dsObs
