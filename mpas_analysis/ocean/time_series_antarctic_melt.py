import os
import xarray

# SFP: the following are only needed for the stop-gap local plotting routine used here
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon

from ..shared.analysis_task import AnalysisTask

from ..shared.constants import constants

from ..shared.plot.plotting import timeseries_analysis_plot

from ..shared.plot.plotting import _plot_xtick_format

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.io.utility import build_config_full_path, make_directories

import csv	

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


        # Load observations from multiple files and put in dictionary based on shelf keyname
        observationsDirectory = build_config_full_path(config, 'oceanObservations', 'meltSubdirectory')
        obsFileNameDict = { 'Rignot et al. (2013)': 'Rignot_2013_melt_rates.csv', 'Rignot et al. (2013) SS': 'Rignot_2013_melt_rates_SS.csv' } 
        obsFileNameList = obsFileNameDict.values() 
        obsFileNameDesc = obsFileNameDict.keys() 

        obsDictDict = {}	# dict for storing dict of obs data
        fileCount = 0		# counter for uniquely identifying obs data set
        while len(obsFileNameList)>0:
	    obsFileName = '{}/{}'.format( observationsDirectory, obsFileNameList.pop(0) )
            fileCount = fileCount + 1
            obsDictTemp = {}
            obsFile = csv.reader( open(obsFileName, 'rU') )	 
            next(obsFile, None)  # skip the header line
    	    for line in obsFile:		# some possibly useful values are left commented out for now
                shelfName = line[0]
                #surveyArea = line[1]
                meltFlux = float( line[2] )
                meltFluxUncertainty = float( line[3] )
                meltRate = float( line[4] )
                meltRateUncertainty = float( line[5] )
                #actualArea = float( line[6] )  # actual area here is in sq km

                # build dict of obs. keyed to filename description (which will be used for plotting)
                obsDictTemp[shelfName] = { 'meltFlux': meltFlux, 'meltFluxUncertainty': meltFluxUncertainty, 
		    'meltRate': meltRate, 'meltRateUncertainty': meltRateUncertainty }
                obsDictDict[ '{}'.format( obsFileNameDesc[fileCount-1] ) ] = obsDictTemp        

	# Note that if areas from obs file are used they need to be converted from sq km to sq m

        # work on data from simulations
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


        obsCount = len( obsDictDict )  
        print '  Make plots...'
        for iRegion in range(nRegions):

            regionName = iceShelvesToPlot[iRegion]

            # get obs melt flux and obs melt flux unc. for shelf (similar for rates) 
            obsMeltFlux = []; obsMeltFluxUnc = []; obsMeltRate = []; obsMeltRateUnc = []
            for iObs in range( obsCount ):
                dictName = '{}'.format( obsFileNameDesc[iObs] )
                obsMeltFlux.append( obsDictDict[dictName][regionName]['meltFlux'] )
                obsMeltFluxUnc.append( obsDictDict[dictName][regionName]['meltFluxUncertainty'] )
                obsMeltRate.append( obsDictDict[dictName][regionName]['meltRate'] )
                obsMeltRateUnc.append( obsDictDict[dictName][regionName]['meltRateUncertainty'] )

            title = regionName.replace('_', ' ')

            regionName = regionName.replace(' ', '_')

            xLabel = 'Time (yr)'
            yLabel = 'Melt Flux (GT/yr)'

            timeSeries = totalMeltFlux.isel(nRegions=iRegion)


            figureName = '{}/melt_flux_{}.png'.format(self.plotsDirectory,
                                                         regionName)

#            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
#                                     title, xLabel, yLabel, figureName,
#                                     lineStyles=['b-'], lineWidths=[1.2],
#                                     calendar=calendar)
            self.plot(config, timeSeries, obsFileNameDesc,  obsCount, obsMeltFlux, obsMeltFluxUnc, title, xLabel, yLabel, figureName )
            


            xLabel = 'Time (yr)'
            yLabel = 'Melt Rate (m/yr)'

            timeSeries = meltRates.isel(nRegions=iRegion)

            figureName = '{}/melt_rate_{}.png'.format(self.plotsDirectory,
                                                         regionName)

#            timeseries_analysis_plot(config, [timeSeries], movingAverageMonths,
#                                     title, xLabel, yLabel, figureName,
#                                     lineStyles=['b-'], lineWidths=[1.2],
#                                     calendar=calendar)
            self.plot(config, timeSeries, obsFileNameDesc, obsCount, obsMeltRate, obsMeltRateUnc, title, xLabel, yLabel, figureName )


    def plot(self, config, modelValues, obsFileNameDesc, obsCount, obsMean, obsUncertainty, title, xlabel, ylabel, fileout):

        """
        Plots the list of time series data sets and stores the result in an image
        file.
        Parameters
        ----------
        config : instance of ConfigParser
            the configuration, containing a [plot] section with options that
            control plotting

        modelValues : xarray DataSet
            the data set to be plotted

        obsMean : float
            a single observed mean value to plot as a constant line

        obsUncertainty : float
            The observed uncertainty, to plot as a shaded rectangle around the mean 

        title : str
            the title of the plot

        xlabel, ylabel : str
            axis labels

        fileout : str
            the file name to be written

        Authors
        -------
        Xylar Asay-Davis, Stephen Price
        """

        # play with these as you need
        figsize = (15, 6)
        dpi = 300
        modelLineStyle = 'b-'
        modelLineWidth = 2
        obsColor = 'gray'
        obsLineWidth = 1 

        maxXTicks = 20
        calendar = self.calendar

        axis_font = {'size': config.get('plot', 'axisFontSize')}
        title_font = {'size': config.get('plot', 'titleFontSize'),
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}

        plt.figure(figsize=figsize, dpi=dpi)

        minDays = modelValues.Time.min()
        maxDays = modelValues.Time.max()
        plt.plot(modelValues.Time.values, modelValues.values, modelLineStyle,
                 linewidth=modelLineWidth, label='model')

        ax = plt.gca()

        # this makes a "patch" with a single rectangular polygon, where the pairs are time and melt 
        # rate for the 4 corners, and "True" means it is a closed polygon:
#        patches = [Polygon([[minDays, obsMean-obsUncertainty], [maxDays, obsMean-obsUncertainty],
#                        [maxDays, obsMean+obsUncertainty], [minDays, obsMean+obsUncertainty]])]

        # make the polygon gray and mostly transparent
#        p = PatchCollection(patches, color=obsColor, alpha=0.15)
        # add it to the plot on the current axis
#        ax.add_collection(p)

        # also plot a line
#        plt.plot([minDays, maxDays], [obsMean, obsMean], color=obsColor, linewidth=obsLineWidth)

	# plot error bars rather than the "patch" used above
        symbol = ['o','^','s','D','*']
        for iObs in range( obsCount ):
           plt.errorbar( ((maxDays - minDays)/5)*(iObs+1)+minDays, obsMean[iObs], yerr=obsUncertainty[iObs], fmt=symbol[iObs], ecolor='k', 
                       capthick=2, label='{}'.format( obsFileNameDesc[iObs] ) )
        # add legend
	plt.legend( loc='lower right', numpoints=1 )

        # this will need to be imported from shared.plot.plotting
        _plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

        if title is not None:
            plt.title(title, **title_font)
        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)
        if fileout is not None:
            plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

        if not config.getboolean('plot', 'displayToScreen'):
            plt.close()

        # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
