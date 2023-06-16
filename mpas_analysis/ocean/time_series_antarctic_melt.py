# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
import os
import xarray
import numpy
import csv
import matplotlib.pyplot as plt

from geometric_features import FeatureCollection, read_feature_collection
from geometric_features.aggregation import get_aggregator_by_name

from mpas_analysis.shared.analysis_task import AnalysisTask

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.plot import timeseries_analysis_plot, savefig, \
    add_inset

from mpas_analysis.shared.io import open_mpas_dataset, write_netcdf_with_fill

from mpas_analysis.shared.io.utility import build_config_full_path, \
    make_directories, build_obs_path, decode_strings

from mpas_analysis.shared.html import write_image_xml


class TimeSeriesAntarcticMelt(AnalysisTask):
    """
    Performs analysis of the time-series output of Antarctic sub-ice-shelf
    melt rates.
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, config, mpasTimeSeriesTask, regionMasksTask,
                 controlConfig=None):

        """
        Construct the analysis task.

        Parameters
        ----------
        config : mpas_tools.config.MpasConfigParser
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        regionMasksTask : ``ComputeRegionMasks``
            A task for computing region masks

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesAntarcticMelt, self).__init__(
            config=config,
            taskName='timeSeriesAntarcticMelt',
            componentName='ocean',
            tags=['timeSeries', 'melt', 'landIceCavities', 'antarctic'])

        regionGroup = 'Ice Shelves'
        iceShelvesToPlot = config.getexpression('timeSeriesAntarcticMelt',
                                                'iceShelvesToPlot')
        if len(iceShelvesToPlot) == 0:
            # nothing else to do
            return

        masksSubtask = regionMasksTask.add_mask_subtask(regionGroup=regionGroup)
        self.iceShelfMasksFile = masksSubtask.geojsonFileName

        iceShelvesToPlot = masksSubtask.expand_region_names(iceShelvesToPlot)

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.get('timeSeries', 'endYear')
        if endYear == 'end':
            # a valid end year wasn't found, so likely the run was not found,
            # perhaps because we're just listing analysis tasks
            endYear = startYear
        else:
            endYear = int(endYear)

        years = list(range(startYear, endYear + 1))

        # in the end, we'll combine all the time series into one, but we
        # create this task first so it's easier to tell it to run after all
        # the compute tasks
        combineSubtask = CombineMeltSubtask(
            self, startYears=years, endYears=years)

        # run one subtask per year
        for year in years:
            computeSubtask = ComputeMeltSubtask(
                self, startYear=year, endYear=year,
                mpasTimeSeriesTask=mpasTimeSeriesTask,
                masksSubtask=masksSubtask,
                iceShelvesToPlot=iceShelvesToPlot)
            self.add_subtask(computeSubtask)
            computeSubtask.run_after(masksSubtask)
            combineSubtask.run_after(computeSubtask)

        self.add_subtask(combineSubtask)

        for index, iceShelf in enumerate(iceShelvesToPlot):
            plotMeltSubtask = PlotMeltSubtask(self, iceShelf, index,
                                              controlConfig)
            plotMeltSubtask.run_after(combineSubtask)
            self.add_subtask(plotMeltSubtask)


class ComputeMeltSubtask(AnalysisTask):
    """
    Computes time-series of Antarctic sub-ice-shelf melt rates.

    Attributes
    ----------
    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    masksSubtask : ``ComputeRegionMasksSubtask``
        A task for creating mask files for each ice shelf to plot

    iceShelvesToPlot : list of str
        A list of ice shelves to plot
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, startYear, endYear, mpasTimeSeriesTask,
                 masksSubtask, iceShelvesToPlot):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  TimeSeriesAntarcticMelt
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        masksSubtask : ``ComputeRegionMasksSubtask``
            A task for creating mask files for each ice shelf to plot

        iceShelvesToPlot : list of str
            A list of ice shelves to plot
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(ComputeMeltSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=f'computeMeltRates_{startYear:04d}-{endYear:04d}')

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.run_after(mpasTimeSeriesTask)

        self.masksSubtask = masksSubtask
        self.run_after(masksSubtask)

        self.iceShelvesToPlot = iceShelvesToPlot
        self.restartFileName = None
        self.startYear = startYear
        self.endYear = endYear
        self.startDate = f'{self.startYear:04d}-01-01_00:00:00'
        self.endDate = f'{self.endYear:04d}-12-31_23:59:59'
        self.variableList = \
            ['timeMonthly_avg_landIceFreshwaterFlux']

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If a restart file is not present

        ValueError
            If ``config_land_ice_flux_mode`` is not one of ``standalone`` or
            ``coupled``
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(ComputeMeltSubtask, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['data', 'standalone', 'coupled']:
            raise ValueError('*** timeSeriesAntarcticMelt requires '
                             'config_land_ice_flux_mode \n'
                             '    to be data, standalone or coupled. '
                             '    Otherwise, no melt rates are available \n'
                             '    for plotting.')

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for Antarctic melt calculations')

        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        return

    def run_task(self):
        """
        Computes time-series of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info("Computing Antarctic melt rate time series...")

        mpasTimeSeriesTask = self.mpasTimeSeriesTask
        config = self.config

        timeSeriesBase = build_config_full_path(config, 'output',
                                                'timeseriesSubdirectory')
        outputDirectory = f'{timeSeriesBase}/iceShelfFluxes/'

        try:
            os.makedirs(outputDirectory)
        except OSError:
            pass

        outFileName = f'{outputDirectory}/iceShelfFluxes_' \
                      f'{self.startYear:04d}-{self.endYear:04d}.nc'

        # Load data:
        inputFile = mpasTimeSeriesTask.outputFile
        dsIn = open_mpas_dataset(fileName=inputFile,
                                 calendar=self.calendar,
                                 variableList=self.variableList,
                                 startDate=self.startDate,
                                 endDate=self.endDate)
        try:
            if os.path.exists(outFileName):
                # The file already exists so load it
                dsOut = xarray.open_dataset(outFileName)
                if numpy.all(dsOut.Time.values == dsIn.Time.values):
                    return
                else:
                    self.logger.warning(f'File {outFileName} is incomplete. '
                                        f'Deleting it.')
                    os.remove(outFileName)
        except OSError:
            # something is potentially wrong with the file, so let's delete
            # it and try again
            self.logger.warning(f'Problems reading file {outFileName}. '
                                f'Deleting it.')
            os.remove(outFileName)

        restartFileName = \
            mpasTimeSeriesTask.runStreams.readpath('restart')[0]

        dsRestart = xarray.open_dataset(restartFileName)
        landIceFraction = dsRestart.landIceFraction.isel(Time=0)
        areaCell = dsRestart.areaCell

        regionMaskFileName = self.masksSubtask.maskFileName

        dsRegionMask = xarray.open_dataset(regionMaskFileName)

        # figure out the indices of the regions to plot
        regionNames = decode_strings(dsRegionMask.regionNames)

        regionIndices = []
        for iceShelf in self.iceShelvesToPlot:
            for index, regionName in enumerate(regionNames):
                if iceShelf == regionName:
                    regionIndices.append(index)
                    break

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)

        regionNames = decode_strings(dsRegionMask.regionNames)

        datasets = []
        nTime = dsIn.sizes['Time']
        for tIndex in range(nTime):
            self.logger.info(f'  {tIndex + 1}/{nTime}')

            freshwaterFlux = \
                dsIn.timeMonthly_avg_landIceFreshwaterFlux.isel(Time=tIndex)

            nRegions = dsRegionMask.sizes['nRegions']
            meltRates = numpy.zeros((nRegions,))
            totalMeltFluxes = numpy.zeros((nRegions,))

            for regionIndex in range(nRegions):
                self.logger.info(f'    {regionNames[regionIndex]}')
                cellMask = \
                    dsRegionMask.regionCellMasks.isel(nRegions=regionIndex)

                # convert from kg/s to kg/yr
                totalMeltFlux = constants.sec_per_year * \
                    (cellMask * areaCell * freshwaterFlux).sum(dim='nCells')

                totalArea = \
                    (landIceFraction * cellMask * areaCell).sum(dim='nCells')

                # from kg/m^2/yr to m/yr
                meltRates[regionIndex] = ((1. / constants.rho_fw) *
                                          (totalMeltFlux / totalArea))

                # convert from kg/yr to GT/yr
                totalMeltFlux /= constants.kg_per_GT
                totalMeltFluxes[regionIndex] = totalMeltFlux

            dsOut = xarray.Dataset()
            dsOut.coords['Time'] = dsIn.Time.isel(Time=tIndex)
            dsOut['totalMeltFlux'] = (('nRegions',), totalMeltFluxes)
            dsOut['meltRates'] = (('nRegions',), meltRates)
            datasets.append(dsOut)

        dsOut = xarray.concat(objs=datasets, dim='Time')
        dsOut['regionNames'] = dsRegionMask.regionNames
        dsOut.totalMeltFlux.attrs['units'] = 'GT a$^{-1}$'
        dsOut.totalMeltFlux.attrs['description'] = \
            'Total melt flux summed over each ice shelf or region'
        dsOut.meltRates.attrs['units'] = 'm a$^{-1}$'
        dsOut.meltRates.attrs['description'] = \
            'Melt rate averaged over each ice shelf or region'

        write_netcdf_with_fill(dsOut, outFileName)


class CombineMeltSubtask(AnalysisTask):
    """
    Combine individual time series into a single data set
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, parentTask, startYears, endYears):
        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask : TimeSeriesAntarcticMelt
            The main task of which this is a subtask

        startYears, endYears : list
            The beginning and end of each time series to combine
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        subtaskName = 'combineAntarcticMeltTimeSeries'

        # first, call the constructor from the base class (AnalysisTask)
        super(CombineMeltSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=subtaskName)

        self.startYears = startYears
        self.endYears = endYears

    def run_task(self):
        """
        Combine the time series
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        timeSeriesBase = build_config_full_path(self.config, 'output',
                                                'timeseriesSubdirectory')
        outputDirectory = f'{timeSeriesBase}/iceShelfFluxes/'

        outFileName = f'{outputDirectory}/iceShelfFluxes_' \
                      f'{self.startYears[0]:04d}-{self.endYears[-1]:04d}.nc'

        if not os.path.exists(outFileName):
            inFileNames = []
            for startYear, endYear in zip(self.startYears, self.endYears):
                inFileName = f'{outputDirectory}/iceShelfFluxes_' \
                             f'{startYear:04d}-{endYear:04d}.nc'
                inFileNames.append(inFileName)

            ds = xarray.open_mfdataset(inFileNames, combine='nested',
                                       concat_dim='Time', decode_times=False)

            ds.load()

            write_netcdf_with_fill(ds, outFileName)


class PlotMeltSubtask(AnalysisTask):
    """
    Plots time-series output of Antarctic sub-ice-shelf melt rates.

    Attributes
    ----------
    iceShelf : str
        Name of the ice shelf to plot

    regionIndex : int
        The index into the dimension ``nRegions`` of the ice shelf to plot

    controlConfig : mpas_tools.config.MpasConfigParser
        The configuration options for the control run (if any)

    """
    # Authors
    # -------
    # Xylar Asay-Davis, Stephen Price

    def __init__(self, parentTask, iceShelf, regionIndex, controlConfig):

        """
        Construct the analysis task.

        Parameters
        ----------
        parentTask :  TimeSeriesAntarcticMelt
            The parent task, used to get the ``taskName``, ``config`` and
            ``componentName``

        iceShelf : str
            Name of the ice shelf to plot

        regionIndex : int
            The index into the dimension ``nRegions`` of the ice shelf to plot

        controlConfig : mpas_tools.config.MpasConfigParser, optional
            Configuration options for a control run (if any)
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(PlotMeltSubtask, self).__init__(
            config=parentTask.config,
            taskName=parentTask.taskName,
            componentName=parentTask.componentName,
            tags=parentTask.tags,
            subtaskName=f'plotMeltRates_{iceShelf.replace(" ", "_")}')

        self.iceShelfMasksFile = parentTask.iceShelfMasksFile
        self.iceShelf = iceShelf
        self.regionIndex = regionIndex
        self.controlConfig = controlConfig

    def setup_and_check(self):
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(PlotMeltSubtask, self).setup_and_check()

        self.xmlFileNames = []

        for prefix in ['melt_flux', 'melt_rate']:
            iceShelfSuffix = self.iceShelf.replace(" ", "_")
            self.xmlFileNames.append(
                f'{self.plotsDirectory}/{prefix}_{iceShelfSuffix}.xml')
        return

    def run_task(self):
        """
        Plots time-series output of Antarctic sub-ice-shelf melt rates.
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Stephen Price

        self.logger.info(f'\nPlotting Antarctic melt rate time series for '
                         f'{self.iceShelf}...')

        self.logger.info('  Load melt rate data...')

        config = self.config
        calendar = self.calendar

        iceShelfMasksFile = self.iceShelfMasksFile

        fcAll = read_feature_collection(iceShelfMasksFile)

        fc = FeatureCollection()
        for feature in fcAll.features:
            if feature['properties']['name'] == self.iceShelf:
                fc.add_feature(feature)
                break

        totalMeltFlux, meltRates = self._load_ice_shelf_fluxes(config)

        plotControl = self.controlConfig is not None
        if plotControl:
            controlRunName = self.controlConfig.get('runs', 'mainRunName')

            refTotalMeltFlux, refMeltRates = \
                self._load_ice_shelf_fluxes(self.controlConfig)
        else:
            controlRunName = None
            refTotalMeltFlux = None
            refMeltRates = None

        # Load observations from multiple files and put in dictionary based
        # on shelf key name
        observationsDirectory = build_obs_path(config, 'ocean',
                                               'meltSubdirectory')
        obsFileNameDict = {'Rignot et al. (2013)':
                           'Rignot_2013_melt_rates_20201117.csv',
                           'Rignot et al. (2013) SS':
                           'Rignot_2013_melt_rates_SS_20201117.csv'}

        obsDict = {}  # dict for storing dict of obs data
        for obsName in obsFileNameDict:
            obsFileName = f'{observationsDirectory}/{obsFileNameDict[obsName]}'
            obsDict[obsName] = {}
            obsFile = csv.reader(open(obsFileName, 'r'))
            next(obsFile, None)  # skip the header line
            for line in obsFile:  # some later useful values commented out
                shelfName = line[0]
                if shelfName != self.iceShelf:
                    continue

                # surveyArea = line[1]
                meltFlux = float(line[2])
                meltFluxUncertainty = float(line[3])
                meltRate = float(line[4])
                meltRateUncertainty = float(line[5])
                # actualArea = float( line[6] )  # actual area here is in sq km

                # build dict of obs. keyed to filename description
                # (which will be used for plotting)
                obsDict[obsName] = {
                    'meltFlux': meltFlux,
                    'meltFluxUncertainty': meltFluxUncertainty,
                    'meltRate': meltRate,
                    'meltRateUncertainty': meltRateUncertainty}
                break
        regionGroup = 'Ice Shelves'
        _, prefix, date = get_aggregator_by_name(regionGroup)

        obsFileName = f'{observationsDirectory}/Adusumilli/Adusumilli_2020_' \
                      f'iceshelf_melt_rates_2010-2018_v0.20230504.' \
                      f'{prefix}{date}.nc'
        with xarray.open_dataset(obsFileName) as ds_adusumilli:
            region_names = [name.values for name in ds_adusumilli.regionNames]
            index = region_names.index(self.iceShelf)
            ds_shelf = ds_adusumilli.isel(nRegions=index)
            obsDict['Adusumilli et al. (2020)'] = {
                'meltFlux': ds_shelf.totalMeltFlux.values,
                'meltFluxUncertainty': ds_shelf.meltFluxUncertainty.values,
                'meltRate': ds_shelf.meanMeltRate.values,
                'meltRateUncertainty': ds_shelf.meltRateUncertainty.values}

        mainRunName = config.get('runs', 'mainRunName')
        movingAveragePoints = config.getint('timeSeriesAntarcticMelt',
                                            'movingAveragePoints')

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        self.logger.info('  Make plots...')

        # get obs melt flux and unc. for shelf (similar for rates)
        obsMeltFlux = []
        obsMeltFluxUnc = []
        obsMeltRate = []
        obsMeltRateUnc = []
        for obsName in obsDict:
            if len(obsDict[obsName]) > 0:
                obsMeltFlux.append(
                    obsDict[obsName]['meltFlux'])
                obsMeltFluxUnc.append(
                    obsDict[obsName]['meltFluxUncertainty'])
                obsMeltRate.append(
                    obsDict[obsName]['meltRate'])
                obsMeltRateUnc.append(
                    obsDict[obsName]['meltRateUncertainty'])
            else:
                # append NaN so this particular obs won't plot
                self.logger.warning(f'{obsName} observations not available '
                                    f'for {self.iceShelf}')
                obsMeltFlux.append(None)
                obsMeltFluxUnc.append(None)
                obsMeltRate.append(None)
                obsMeltRateUnc.append(None)

        title = self.iceShelf.replace('_', ' ')
        suffix = self.iceShelf.replace(' ', '_')

        xLabel = 'Time (yr)'
        yLabel = 'Melt Flux (GT/yr)'

        timeSeries = totalMeltFlux.isel(nRegions=self.regionIndex)

        filePrefix = f'melt_flux_{suffix}'
        outFileName = f'{self.plotsDirectory}/{filePrefix}.png'

        fields = [timeSeries]
        lineColors = [config.get('timeSeries', 'mainColor')]
        lineWidths = [2.5]
        legendText = [mainRunName]
        if plotControl:
            fields.append(refTotalMeltFlux.isel(nRegions=self.regionIndex))
            lineColors.append(config.get('timeSeries', 'controlColor'))
            lineWidths.append(1.2)
            legendText.append(controlRunName)

        if config.has_option('timeSeriesAntarcticMelt', 'firstYearXTicks'):
            firstYearXTicks = config.getint('timeSeriesAntarcticMelt',
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option('timeSeriesAntarcticMelt', 'yearStrideXTicks'):
            yearStrideXTicks = config.getint('timeSeriesAntarcticMelt',
                                             'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        if config.has_option('timeSeriesAntarcticMelt', 'titleFontSize'):
            titleFontSize = config.getint('timeSeriesAntarcticMelt',
                                          'titleFontSize')
        else:
            titleFontSize = None

        if config.has_option('timeSeriesAntarcticMelt', 'defaultFontSize'):
            defaultFontSize = config.getint('timeSeriesAntarcticMelt',
                                            'defaultFontSize')
        else:
            defaultFontSize = None

        fig = timeseries_analysis_plot(config, fields, calendar=calendar,
                                       title=title, xlabel=xLabel,
                                       ylabel=yLabel,
                                       movingAveragePoints=movingAveragePoints,
                                       lineColors=lineColors,
                                       lineWidths=lineWidths,
                                       legendText=legendText,
                                       legendLocation='upper left',
                                       titleFontSize=titleFontSize,
                                       defaultFontSize=defaultFontSize,
                                       obsMean=obsMeltFlux,
                                       obsUncertainty=obsMeltFluxUnc,
                                       obsLegend=list(obsDict.keys()),
                                       firstYearXTicks=firstYearXTicks,
                                       yearStrideXTicks=yearStrideXTicks)

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName, config)

        caption = f'Running Mean of Total Melt Flux under Ice Shelves in ' \
                  f'the {title} Region'
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Antarctic Melt Time Series',
            groupLink='antmelttime',
            gallery='Total Melt Flux',
            thumbnailDescription=title,
            imageDescription=caption,
            imageCaption=caption)

        xLabel = 'Time (yr)'
        yLabel = 'Melt Rate (m/yr)'

        timeSeries = meltRates.isel(nRegions=self.regionIndex)

        filePrefix = f'melt_rate_{suffix}'
        outFileName = f'{self.plotsDirectory}/{filePrefix}.png'

        fields = [timeSeries]
        lineColors = [config.get('timeSeries', 'mainColor')]
        lineWidths = [2.5]
        legendText = [mainRunName]
        if plotControl:
            fields.append(refMeltRates.isel(nRegions=self.regionIndex))
            lineColors.append(config.get('timeSeries', 'controlColor'))
            lineWidths.append(1.2)
            legendText.append(controlRunName)

        if config.has_option(self.taskName, 'firstYearXTicks'):
            firstYearXTicks = config.getint(self.taskName,
                                            'firstYearXTicks')
        else:
            firstYearXTicks = None

        if config.has_option(self.taskName, 'yearStrideXTicks'):
            yearStrideXTicks = config.getint(self.taskName,
                                             'yearStrideXTicks')
        else:
            yearStrideXTicks = None

        fig = timeseries_analysis_plot(config, fields, calendar=calendar,
                                       title=title, xlabel=xLabel,
                                       ylabel=yLabel,
                                       movingAveragePoints=movingAveragePoints,
                                       lineColors=lineColors,
                                       lineWidths=lineWidths,
                                       legendText=legendText,
                                       firstYearXTicks=firstYearXTicks,
                                       yearStrideXTicks=yearStrideXTicks,
                                       obsMean=obsMeltRate,
                                       obsUncertainty=obsMeltRateUnc,
                                       obsLegend=list(obsDict.keys()))

        # do this before the inset because otherwise it moves the inset
        # and cartopy doesn't play too well with tight_layout anyway
        plt.tight_layout()

        add_inset(fig, fc, width=2.0, height=2.0)

        savefig(outFileName, config)

        caption = f'Running Mean of Area-averaged Melt Rate under Ice ' \
                  f'Shelves in the {title} Region'
        write_image_xml(
            config=config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup='Antarctic Melt Time Series',
            groupLink='antmelttime',
            gallery='Area-averaged Melt Rate',
            thumbnailDescription=title,
            imageDescription=caption,
            imageCaption=caption)

    @staticmethod
    def _load_ice_shelf_fluxes(config):
        """
        Reads melt flux time series and computes regional total melt flux and
        mean melt rate.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        timeSeriesBase = build_config_full_path(config, 'output',
                                                'timeseriesSubdirectory')

        outputDirectory = f'{timeSeriesBase}/iceShelfFluxes/'

        startYear = config.getint('timeSeries', 'startYear')
        endYear = config.getint('timeSeries', 'endYear')

        outFileName = f'{outputDirectory}/iceShelfFluxes_' \
                      f'{startYear:04d}-{endYear:04d}.nc'

        dsOut = xarray.open_dataset(outFileName)
        return dsOut.totalMeltFlux, dsOut.meltRates
