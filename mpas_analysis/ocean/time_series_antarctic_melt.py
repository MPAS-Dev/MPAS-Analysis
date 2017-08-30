import os
import xarray

from ..shared.analysis_task import AnalysisTask

from ..shared.constants import constants

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.io.utility import build_config_full_path, make_directories


class TimeSeriesAntarcticMelt(AnalysisTask):
    """
    Performs analysis of the time-series output of Antarctic sub-ice-shelf
    melt rates.

    Authors
    -------
    Xylar Asay-Davis, Stephen Price
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
        # first, call the constructor from the base class (AnalysisTask)
        super(TimeSeriesAntarcticMelt, self).__init__(
            config=config,
            taskName='timeSeriesAntarcticMelt',
            componentName='ocean',
            tags=['timeSeries', 'melt', 'landIceCavities'])

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
        streamName = self.historyStreams.find_stream(
            self.streamMap['timeSeriesStats'])
        self.startDate = config.get('timeSeries', 'startDate')
        self.endDate = config.get('timeSeries', 'endDate')
        self.inputFiles = \
            self.historyStreams.readpath(streamName,
                                         startDate=self.startDate,
                                         endDate=self.endDate,
                                         calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDate,
                                       self.endDate))

        return  # }}}

    def run(self):  # {{{
        """
        Performs analysis of the time-series output of Antarctic sub-ice-shelf
        melt rates.

        Authors
        -------
        Xylar Asay-Davis
        """

        print "\nPlotting Antarctic melt rate time series..."

        print '  Load melt rate data...'

        config = self.config
        calendar = self.calendar

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))

        # Load data:
        variableList = ['timeMonthly_avg_landIceFreshwaterFlux']
        timeVariableName = ['xtime_startMonthly', 'xtime_endMonthly']
        ds = open_multifile_dataset(fileNames=self.inputFiles,
                                    calendar=calendar,
                                    config=config,
                                    timeVariableName=timeVariableName,
                                    variableList=variableList,
                                    startDate=self.startDate,
                                    endDate=self.endDate)

        freshwaterFlux = ds.timeMonthly_avg_landIceFreshwaterFlux

        movingAverageMonths = config.getint('timeSeriesAntarcticMelt',
                                            'movingAverageMonths')

        iceShelvesToPlot = config.getExpression('timeSeriesAntarcticMelt',
                                                'iceShelvesToPlot')

        dsRestart = xarray.open_dataset(self.restartFileName)
        areaCell = dsRestart.landIceFraction.isel(Time=0)*dsRestart.areaCell

        dsRegionMask = xarray.open_dataset(self.regionMaskFileName)
        regionNames = list(dsRegionMask.regionNames.values)
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

        # select only those regions we want to plot
        dsRegionMask = dsRegionMask.isel(nRegions=regionIndices)
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

        print '  Make plots...'
        for iRegion in range(nRegions):
            regionName = iceShelvesToPlot[iRegion]

            title = regionName.replace('_', ' ')

            regionName = regionName.replace(' ', '_')

            xLabel = 'Time (yr)'
            yLabel = 'Melt Flux (GT/yr)'

            timeSeries = totalMeltFlux.isel(nRegions=iRegion)

            figureName = '{}/melt_flux_{}.png'.format(self.plotsDirectory,
                                                         regionName)

            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     lineStyles=['b-'], lineWidths=[1.2],
                                     calendar=calendar)

            xLabel = 'Time (yr)'
            yLabel = 'Melt Rate (m/yr)'

            timeSeries = meltRates.isel(nRegions=iRegion)

            figureName = '{}/melt_rate_{}.png'.format(self.plotsDirectory,
                                                         regionName)

            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
                                     title, xLabel, yLabel, figureName,
                                     lineStyles=['b-'], lineWidths=[1.2],
                                     calendar=calendar)
        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
