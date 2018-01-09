from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import xarray


from ..shared.analysis_task import AnalysisTask

from ..shared.constants import constants

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.io import open_mpas_dataset

from ..shared.io.utility import build_config_full_path, make_directories

from ..shared.html import write_image_xml

import csv


class TimeSeriesAntarcticMelt(AnalysisTask):
    """
    Performs analysis of the time-series output of Antarctic sub-ice-shelf
    melt rates.

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    Authors
    -------
    Xylar Asay-Davis, Stephen Price
    """

    def __init__(self, config, mpasTimeSeriesTask):  # {{{
        """
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesAntarcticMelt, self).__init__(
            config=config,
            taskName='timeSeriesAntarcticMelt',
            componentName='ocean',
            tags=['timeSeries', 'melt', 'landIceCavities'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask

        self.run_after(mpasTimeSeriesTask)

        # }}}

    def setup_and_check(self):  # {{{
        """
        Perform steps to set up the analysis and check for errors in the setup.

        Raises
        ------
        IOError
            If files are not present

        Authors
        -------
        Xylar Asay-Davis
        """
        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #   self.inDirectory, self.plotsDirectory, self.namelist, self.streams
        #   self.calendar
        super(TimeSeriesAntarcticMelt, self).setup_and_check()

        self.check_analysis_enabled(
            analysisOptionName='config_am_timeseriesstatsmonthly_enable',
            raiseException=True)

        config = self.config

        landIceFluxMode = self.namelist.get('config_land_ice_flux_mode')
        if landIceFluxMode not in ['standalone', 'coupled']:
            raise ValueError('*** timeSeriesAntarcticMelt requires '
                             'config_land_ice_flux_mode \n'
                             '    to be standalone or coupled.  Otherwise, no '
                             'melt rates are available \n'
                             '    for plotting.')

        mpasMeshName = config.get('input', 'mpasMeshName')
        regionMaskDirectory = config.get('regions', 'regionMaskDirectory')

        self.regionMaskFileName = '{}/{}_iceShelfMasks.nc'.format(
                regionMaskDirectory, mpasMeshName)

        if not os.path.exists(self.regionMaskFileName):
            raise IOError('Regional masking file {} for Antarctica melt-rate '
                          'calculation does not exist'.format(
                                  self.regionMaskFileName))

        # Load mesh related variables
        try:
            self.restartFileName = self.runStreams.readpath('restart')[0]
        except ValueError:
            raise IOError('No MPAS-O restart file found: need at least one '
                          'restart file for Antarctic melt calculations')

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')

        self.inputFile = self.mpasTimeSeriesTask.outputFile

        self.variableList = \
            ['timeMonthly_avg_landIceFreshwaterFlux']
        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        iceShelvesToPlot = config.getExpression('timeSeriesAntarcticMelt',
                                                'iceShelvesToPlot')

        with xarray.open_dataset(self.regionMaskFileName) as dsRegionMask:
            regionNames = [bytes.decode(name) for name in
                           dsRegionMask.regionNames.values]
            nRegions = dsRegionMask.dims['nRegions']

        if 'all' in iceShelvesToPlot:
            iceShelvesToPlot = regionNames
            regionIndices = [iRegion for iRegion in range(nRegions)]

        else:
            regionIndices = []
            for regionName in iceShelvesToPlot:
                if regionName not in regionNames:
                    raise ValueError('Unknown ice shelf name {}'.format(
                            regionName))

                iRegion = regionNames.index(regionName)
                regionIndices.append(iRegion)

        self.regionIndices = regionIndices
        self.iceShelvesToPlot = iceShelvesToPlot
        self.xmlFileNames = []

        for prefix in ['melt_flux', 'melt_rate']:
            for regionName in iceShelvesToPlot:
                regionName = regionName.replace(' ', '_')
                self.xmlFileNames.append(
                    '{}/{}_{}.xml'.format(self.plotsDirectory, prefix,
                                          regionName))
        return  # }}}

    def run_task(self):  # {{{
        """
        Performs analysis of the time-series output of Antarctic sub-ice-shelf
        melt rates.

        Authors
        -------
        Xylar Asay-Davis
        """

        self.logger.info("\nPlotting Antarctic melt rate time series...")

        self.logger.info('  Load melt rate data...')

        config = self.config
        calendar = self.calendar

        # Load data:
        ds = open_mpas_dataset(fileName=self.inputFile,
                               calendar=calendar,
                               variableList=self.variableList,
                               startDate=self.startDate,
                               endDate=self.endDate)

        # Load observations from multiple files and put in dictionary based
        # on shelf keyname
        observationsDirectory = build_config_full_path(config,
                                                       'oceanObservations',
                                                       'meltSubdirectory')
        obsFileNameDict = {'Rignot et al. (2013)':
                           'Rignot_2013_melt_rates.csv',
                           'Rignot et al. (2013) SS':
                           'Rignot_2013_melt_rates_SS.csv'}

        obsDict = {}  # dict for storing dict of obs data
        for obsName in obsFileNameDict:
            obsFileName = '{}/{}'.format(observationsDirectory,
                                         obsFileNameDict[obsName])
            obsDict[obsName] = {}
            obsFile = csv.reader(open(obsFileName, 'rU'))
            next(obsFile, None)  # skip the header line
            for line in obsFile:  # some later useful values commented out
                shelfName = line[0]
                # surveyArea = line[1]
                meltFlux = float(line[2])
                meltFluxUncertainty = float(line[3])
                meltRate = float(line[4])
                meltRateUncertainty = float(line[5])
                # actualArea = float( line[6] )  # actual area here is in sq km

                # build dict of obs. keyed to filename description
                # (which will be used for plotting)
                obsDict[obsName][shelfName] = {
                        'meltFlux': meltFlux,
                        'meltFluxUncertainty': meltFluxUncertainty,
                        'meltRate': meltRate,
                        'meltRateUncertainty': meltRateUncertainty}

        # If areas from obs file used need to be converted from sq km to sq m

        # work on data from simulations
        freshwaterFlux = ds.timeMonthly_avg_landIceFreshwaterFlux

        mainRunName = config.get('runs', 'mainRunName')
        movingAverageMonths = config.getint('timeSeriesAntarcticMelt',
                                            'movingAverageMonths')

        iceShelvesToPlot = self.iceShelvesToPlot

        dsRestart = xarray.open_dataset(self.restartFileName)
        areaCell = dsRestart.landIceFraction.isel(Time=0)*dsRestart.areaCell

        dsRegionMask = xarray.open_dataset(self.regionMaskFileName)

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=self.regionIndices)
        cellMasks = dsRegionMask.regionCellMasks
        nRegions = dsRegionMask.dims['nRegions']

        # convert from kg/s to kg/yr
        totalMeltFlux = constants.sec_per_year * \
            (cellMasks*areaCell*freshwaterFlux).sum(dim='nCells')

        totalArea = (cellMasks*areaCell).sum(dim='nCells')

        # from kg/m^2/yr to m/yr
        meltRates = (1./constants.rho_fw) * (totalMeltFlux/totalArea)

        # convert from kg/yr to GT/yr
        totalMeltFlux /= constants.kg_per_GT

        outputDirectory = build_config_full_path(config, 'output',
                                                 'timeseriesSubdirectory')

        make_directories(outputDirectory)

        self.logger.info('  Make plots...')
        for iRegion in range(nRegions):

            regionName = iceShelvesToPlot[iRegion]

            # get obs melt flux and unc. for shelf (similar for rates)
            obsMeltFlux = []
            obsMeltFluxUnc = []
            obsMeltRate = []
            obsMeltRateUnc = []
            for obsName in obsDict:
                obsMeltFlux.append(
                    obsDict[obsName][regionName]['meltFlux'])
                obsMeltFluxUnc.append(
                    obsDict[obsName][regionName]['meltFluxUncertainty'])
                obsMeltRate.append(
                    obsDict[obsName][regionName]['meltRate'])
                obsMeltRateUnc.append(
                    obsDict[obsName][regionName]['meltRateUncertainty'])

            title = regionName.replace('_', ' ')

            regionName = regionName.replace(' ', '_')

            xLabel = 'Time (yr)'
            yLabel = 'Melt Flux (GT/yr)'

            timeSeries = totalMeltFlux.isel(nRegions=iRegion)

            filePrefix = 'melt_flux_{}'.format(regionName)
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     lineStyles=['b-'], lineWidths=[1.2],
                                     legendText=[mainRunName],
                                     calendar=calendar, obsMean=obsMeltFlux,
                                     obsUncertainty=obsMeltFluxUnc,
                                     obsLegend=list(obsDict.keys()))

            caption = 'Running Mean of Total Melt Flux  under Ice ' \
                      'Shelves in the {} Region'.format(title)
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

            timeSeries = meltRates.isel(nRegions=iRegion)

            filePrefix = 'melt_rate_{}'.format(regionName)
            figureName = '{}/{}.png'.format(self.plotsDirectory, filePrefix)

            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     lineStyles=['b-'], lineWidths=[1.2],
                                     legendText=[mainRunName],
                                     calendar=calendar, obsMean=obsMeltRate,
                                     obsUncertainty=obsMeltRateUnc,
                                     obsLegend=list(obsDict.keys()))

            caption = 'Running Mean of Area-averaged Melt Rate under Ice ' \
                      'Shelves in the {} Region'.format(title)
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
        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
