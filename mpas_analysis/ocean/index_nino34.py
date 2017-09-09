# -*- coding: utf-8 -*-
import datetime
import xarray as xr
import pandas as pd
import numpy as np
from scipy import signal, stats
import os
import matplotlib.pyplot as plt

from ..shared.climatology import climatology
from ..shared.constants import constants
from ..shared.io.utility import build_config_full_path
from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.plot.plotting import plot_xtick_format, plot_size_y_axis

from ..shared.analysis_task import AnalysisTask
from ..shared.html import write_image_xml


class IndexNino34(AnalysisTask):  # {{{
    '''
    <Describe the analysis task here.>

    Authors
    -------
    <List of authors>
    '''

    def __init__(self, config):  # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  instance of MpasAnalysisConfigParser
            Contains configuration options

        Authors
        -------
        Xylar Asay-Davis

        '''
        # first, call the constructor from the base class (AnalysisTask)
        super(IndexNino34, self).__init__(
            config=config,
            taskName='indexNino34',
            componentName='ocean',
            tags=['index', 'nino'])

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.

        Authors
        -------
        Xylar Asay-Davis
        '''

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(IndexNino34, self).setup_and_check()

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        streamName = 'timeSeriesStatsMonthlyOutput'
        self.startDate = self.config.get('index', 'startDate')
        self.endDate = self.config.get('index', 'endDate')
        self.inputFiles = self.historyStreams.readpath(
                streamName, startDate=self.startDate, endDate=self.endDate,
                calendar=self.calendar)

        if len(self.inputFiles) == 0:
            raise IOError('No files were found in stream {} between {} and '
                          '{}.'.format(streamName, self.startDate,
                                       self.endDate))

        mainRunName = self.config.get('runs', 'mainRunName')

        self.xmlFileNames = []
        for filePrefix in ['NINO34_{}'.format(mainRunName),
                           'NINO34_spectra_{}'.format(mainRunName)]:
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))

        # }}}

    def run(self):  # {{{
        '''
        Computes NINO34 index and plots the time series and power spectrum with
        95 and 99% confidence bounds

        Authors
        -------
        Luke Van Roekel, Xylar Asay-Davis
        '''

        print "\nPlotting Nino3.4 time series and power spectrum...."

        print '  Load SST data...'
        fieldName = 'nino'

        simulationStartTime = get_simulation_start_time(self.runStreams)
        config = self.config
        calendar = self.calendar

        dataSource = config.get('indexNino34', 'observationData')

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', '{}Subdirectory'.format(fieldName))

        # specify obsTitle based on data path
        # These are the only data sets supported
        if dataSource == 'HADIsst':
            dataPath = "{}/HADIsst_nino34.nc".format(observationsDirectory)
            obsTitle = 'HADSST'
        else:
            dataPath = "{}/ERS_SSTv4_nino34.nc".format(observationsDirectory)
            obsTitle = 'ERS SSTv4'

        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(self.inputFiles[0]),
                  os.path.basename(self.inputFiles[-1]))
        mainRunName = config.get('runs', 'mainRunName')

        # regionIndex should correspond to NINO34 in surface weighted Average
        # AM
        regionIndex = config.getint('indexNino34', 'regionIndicesToPlot')

        # Load data:
        varName = \
            'timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature'
        varList = [varName]
        ds = open_multifile_dataset(fileNames=self.inputFiles,
                                    calendar=calendar,
                                    config=config,
                                    simulationStartTime=simulationStartTime,
                                    timeVariableName=['xtime_startMonthly',
                                                      'xtime_endMonthly'],
                                    variableList=varList,
                                    startDate=self.startDate,
                                    endDate=self.endDate)

        # Observations have been processed to the nino34Index prior to reading
        dsObs = xr.open_dataset(dataPath)
        nino34Obs = dsObs.sst

        print '  Compute NINO3.4 index...'
        regionSST = ds[varName].isel(nOceanRegions=regionIndex)
        nino34 = self._compute_nino34_index(regionSST, calendar)

        # Compute the observational index over the entire time range
        # nino34Obs = compute_nino34_index(dsObs.sst, calendar)

        print ' Computing NINO3.4 power spectra...'
        f, spectra, conf99, conf95, redNoise = \
            self._compute_nino34_spectra(nino34)

        # Compute the observational spectra over the whole record
        fObs, spectraObs, conf99Obs, conf95Obs, redNoiseObs = \
            self._compute_nino34_spectra(nino34Obs)

        # Compute the observational spectra over the last 30 years for
        # comparison. Only saving the spectra
        time_start = datetime.datetime(1976, 1, 1)
        time_end = datetime.datetime(2016, 12, 31)
        nino3430 = nino34Obs.sel(Time=slice(time_start, time_end))
        f30, spectra30yrs, conf9930, conf9530, redNoise30 = \
            self._compute_nino34_spectra(nino3430)

        # Convert frequencies to period in years
        f = 1.0 / (constants.eps + f*constants.sec_per_year)
        fObs = 1.0 / (constants.eps + fObs*constants.sec_per_year)
        f30 = 1.0 / (constants.eps + f30*constants.sec_per_year)

        print ' Plot NINO3.4 index and spectra...'

        figureName = '{}/NINO34_{}.png'.format(self.plotsDirectory,
                                               mainRunName)
        modelTitle = "{}".format(mainRunName)
        self._nino34_timeseries_plot(config, nino34, nino34Obs, nino3430,
                                     'NINO 3.4 Index', modelTitle, obsTitle,
                                     figureName, linewidths=2,
                                     calendar=calendar)

        self._write_xml(filePrefix='NINO34_{}'.format(mainRunName),
                        plotType='Time Series')

        figureName = '{}/NINO34_spectra_{}.png'.format(self.plotsDirectory,
                                                       mainRunName)
        self._nino34_spectra_plot(config, f, spectra, conf95, conf99, redNoise,
                                  fObs, f30, spectraObs, conf95Obs, conf99Obs,
                                  redNoiseObs, spectra30yrs, conf9530,
                                  conf9930, redNoise30,
                                  'NINO3.4 power spectrum', modelTitle,
                                  obsTitle, figureName, linewidths=2)

        self._write_xml(filePrefix='NINO34_spectra_{}'.format(mainRunName),
                        plotType='Spectra')

    # }}}

    def _compute_nino34_index(self, regionSST, calendar):  # {{{
        """
        Computes nino34 index time series.  It follow the standard nino34
        algorithm, i.e.,

          1. Compute monthly average SST in the region
          2. Computes anomalous SST
          3. Performs a 5 month running mean over the anomalies

        This routine requires regionSST to be the SSTs in the nino3.4 region
        ONLY. It is defined as lat > -5S and lat < 5N and lon > 190E and
        lon < 240E.

        Parameters
        ----------
        regionSST : xarray.DataArray object
           values of SST in the nino region

        calendar: {'gregorian', 'gregorian_noleap'}
            The name of the calendars used in the MPAS run

        Returns
        -------
        xarray.DataArray object containing the nino34index

        Author
        ------
        Luke Van Roekel, Xylar Asay-Davis
        """

        if not isinstance(regionSST, xr.core.dataarray.DataArray):
            raise ValueError('regionSST should be an xarray DataArray')

        # add 'month' data array so we can group by month below.
        regionSST = climatology.add_years_months_days_in_month(regionSST,
                                                               calendar)

        # Compute monthly average and anomaly of climatology of SST
        monthlyClimatology = \
            climatology.compute_monthly_climatology(regionSST,
                                                    maskVaries=False)

        anomaly = regionSST.groupby('month') - monthlyClimatology

        # Remove the long term trend from the anomalies
        detrendedAnomal = signal.detrend(anomaly.values)
        anomaly.values = detrendedAnomal

        # Compute 5 month running mean
        wgts = np.ones(5) / 5.
        return self._running_mean(anomaly, wgts)  # }}}

    def _compute_nino34_spectra(self, nino34Index):  # {{{
        """
        Computes power spectra of Nino34 index.

        nino34Index is the NINO index computed by compute_nino34_index

        The algorithm follows the NCL cvdp package see
        http://www.cesm.ucar.edu/working_groups/CVC/cvdp/code.html

        Parameters
        ----------
        nino34Index : xarray.DataArray object
            nino34Index for analysis

        Returns
        -------
        pxxSmooth : xarray.DataArray object
            nino34Index power spectra that has been smoothed with a modified
            Daniell window (https://www.ncl.ucar.edu/Document/Functions/Built-in/specx_anal.shtml)


        f : numpy.array
            array of frequencies corresponding to the center of the spectral
            bins resulting from the analysis

        mkov*scale : numpy.array
            Red noise fit to pxxSmooth

        mkov*scale*xLow : numpy.array
            95% confidence threshold from chi-squared test

        mkov*scale*xHigh : numpy.array
            99% confidence threshold from chi-squared test

        Author
        ------
        Luke Van Roekel, Xylar Asay-Davis
        """

        # Move nino34Index to numpy to allow functionality with scipy routines
        ninoIndex = nino34Index.values
        window = signal.tukey(len(ninoIndex), alpha=0.1)
        f, Pxx = signal.periodogram(window * ninoIndex,
                                    1.0 / constants.sec_per_month)

        # computes power spectra, smoothed with a weighted running mean
        nwts = max(1, int(7*len(ninoIndex) / 1200))
        # verify window length is odd, if not, add 1
        if nwts % 2 == 0:
            nwts += 1
        # Calculate the weights for the running mean
        # Weights are from the modified Daniell Window
        wgts = np.ones(nwts)
        wgts[0] = 0.5
        wgts[-1] = 0.5
        wgts /= sum(wgts)

        pxxSmooth = (self._running_mean(pd.Series(Pxx), wgts) /
                     constants.sec_per_month)

        # compute 99 and 95% confidence intervals and red-noise process
        # Uses Chi squared test

        r = self._autocorr(ninoIndex)[0, 1]
        r2 = 2.*r
        rsq = r**2

        # In the temp2 variable, f is converted to give wavenumber, i.e.
        # 0,1,2,...,N/2
        temp2 = r2*np.cos(2.*np.pi*f*constants.sec_per_month)
        mkov = 1. / (1. + rsq - temp2)

        sum1 = np.sum(mkov)
        sum2 = np.sum(pxxSmooth.values)
        scale = sum2 / sum1

        df = 2. / (constants.tapcoef * sum(wgts**2))
        xLow = stats.chi2.interval(0.95, df)[1]/df
        xHigh = stats.chi2.interval(0.99, df)[1]/df

        # return Spectra, 99% confidence level, 95% confidence level,
        #        and Red-noise fit
        return f, pxxSmooth, mkov*scale*xHigh, mkov*scale*xLow, mkov*scale
        # }}}

    def _autocorr(self, x, t=1):  # {{{
        """
        Computes lag one auto-correlation for the NINO34 spectra calculation

        Parameters
        ----------
        x : numpy 1-D array
           time series array

        Returns
        -------
        Single value giving the lag one auto-correlation
            If t != 1, this is no longer a lag one auto-correlation

        Author
        ------
        Luke Van Roekel
        """

        return np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))  # }}}

    def _running_mean(self, inputData, wgts):  # {{{
        """
        Calculates a generic weighted running mean

        Parameters
        ----------
        inputData : xr.DataArray
           Data to be smoothed

        wgts : numpy.array
           array of weights that give the smoothing type
           for the nino index this is a 5-point boxcar window
           for the nino power spectra this is a modified Daniell window (see
           https://www.ncl.ucar.edu/Document/Functions/Built-in/specx_anal.shtml)

        Author
        ------
        Luke Van Roekel, Xylar Asay-Davis
        """

        nt = len(inputData)
        sp = (len(wgts) - 1)/2
        runningMean = inputData.copy()
        for k in range(sp, nt-(sp+1)):
            runningMean[k] = sum(wgts*inputData[k-sp:k+sp+1].values)

        return runningMean  # }}}

    def _nino34_spectra_plot(self, config, f, ninoSpectra,
                             confidence95, confidence99, redNoiseSpectra,
                             fObs, f30, ninoObs,
                             conf95Obs, conf99Obs, redNoiseObs,
                             nino30yr, conf9530, conf9930, redNoise30,
                             title, modelTitle, obsTitle,
                             fileout, linewidths, xlabel='Period (years)',
                             ylabel=r'Power ($^o$C / cycles mo$^{-1}$)',
                             titleFontSize=None, figsize=(9, 21), dpi=None):
        # {{{
        """
        Plots the nino34 time series and power spectra in an image file
        Parameters
        ----------
        config : instance of ConfigParser
            the configuration, containing a [plot] section with options that
            control plotting

        f : numpy.array
            periods to plot on x-axis

        ninoSpectra : xarray.dataArray object
            nino34 power spectra

        confidence95 : numpy.array
            95% confidence level based on chi squared test

        confidence99 : numpy.array
            99% confidence level based on chi squared test

        redNoiseSpectra : numpy.array
            red noise fit to the ninoSpectra

        fObs : numpy.array
               periods to plot on x-axis for observations

        ninoObs : xarray.dataArray object
            nino34 power spectra from the full observational record

        conf95Obs : numpy.array
            95% confidence level based on chi squared for observations

        conf99Obs : numpy.array
            99% confidence level based on chi squared for observations

        redNoiseObs : numpy.array
            red noise fit to ninoObs

        nino30yr : xarray.dataArray object
            power spectra of the last 30 years of the observational record

        title : str
            the title of the plot

        modelTitle : str
            the title of model panel

        obsTitle : str
            the title of the obs panel

        xLabel, yLabel : str
            axis labels

        fileout : str
            the file name to be written

        linewidths : control line width

        titleFontSize : int, optional
            the size of the title font

        figsize : tuple of float, optional
            the size of the figure in inches

        dpi : int, optional
            the number of dots per inch of the figure, taken from section
            ``plot`` option ``dpi`` in the config file by default

        Author
        ------
        Luke Van Roekel, Xylar Asay-Davis
        """

        if dpi is None:
            dpi = config.getint('plot', 'dpi')
        fig = plt.figure(figsize=figsize, dpi=dpi)

        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')

        axis_font = {'size': config.get('plot', 'axisFontSize')}
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        if title is not None:
            fig.suptitle(title, y=0.92, **title_font)

        ax1 = plt.subplot(3, 1, 1)

        plt.plot(fObs[2:-3], ninoObs[2:-3], 'k', linewidth=linewidths)
        plt.plot(fObs[2:-3], redNoiseObs[2:-3], 'r', linewidth=linewidths)
        plt.plot(fObs[2:-3], conf95Obs[2:-3], 'b', linewidth=linewidths)
        plt.plot(fObs[2:-3], conf99Obs[2:-3], 'g', linewidth=linewidths)
        plt.xlim(10, 1)

        plt.legend(['Nino34 spectra (Full Record)', 'Red noise fit',
                   '95% confidence threshold', '99% confidence threshold'],
                   loc='upper right')
        maxObs = plot_size_y_axis(plt, fObs, c1=conf99Obs, c2=redNoiseObs)
        max30 = plot_size_y_axis(plt, f30, c1=conf9930, c2=redNoise30)
        maxModel = plot_size_y_axis(plt, f, c1=ninoSpectra.values,
                                    c2=confidence99, c3=redNoiseSpectra)

        maxYval = max(maxObs, max30, maxModel)
        plt.ylim(0, 0.9*maxYval)

        if obsTitle is not None:
            plt.title(obsTitle+' (Full Record)', **title_font)
        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)

        ax2 = plt.subplot(3, 1, 2)

        plt.plot(f30[2:-3], nino30yr[2:-3], 'k', linewidth=linewidths)
        plt.plot(f30[2:-3], redNoise30[2:-3], 'r', linewidth=linewidths)
        plt.plot(f30[2:-3], conf9530[2:-3], 'b', linewidth=linewidths)
        plt.plot(f30[2:-3], conf9930[2:-3], 'g', linewidth=linewidths)
        plt.xlim(10, 1)
        plt.ylim(0, 0.9*maxYval)

        plt.legend(['Nino34 spectra (1976 - 2016)', 'Red noise fit',
                   '95% confidence threshold', '99% confidence threshold'],
                   loc='upper right')

        if obsTitle is not None:
            plt.title(obsTitle+' (1976-2016)', **title_font)
        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)

        ax3 = plt.subplot(3, 1, 3)
        plt.plot(f[2:-3], ninoSpectra[2:-3], 'k', linewidth=linewidths)
        plt.plot(f[2:-3], redNoiseSpectra[2:-3], 'r', linewidth=linewidths)
        plt.plot(f[2:-3], confidence95[2:-3], 'b', linewidth=linewidths)
        plt.plot(f[2:-3], confidence99[2:-3], 'g', linewidth=linewidths)
        plt.xlim(10, 1)
        plt.ylim(0, 0.9*maxYval)

        # add legend
        plt.legend(['Nino34 index spectra', 'Red noise fit',
                   '95% confidence threshold', '99% confidence threshold'],
                   loc='upper right')

        if modelTitle is not None:
            plt.title(modelTitle, **title_font)
        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)
        if fileout is not None:
            fig.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

        if not config.getboolean('plot', 'displayToScreen'):
            plt.close()
        # }}}

    def _nino34_timeseries_plot(self, config, nino34Index, nino34Obs, nino3430,
                                title, modelTitle, obsTitle, fileout,
                                linewidths, calendar, xlabel='Time [years]',
                                ylabel='[$^\circ$C]', titleFontSize=None,
                                figsize=(12, 28), dpi=None, maxXTicks=20):
        # {{{
        """
        Plots the nino34 time series and power spectra in an image file

        Parameters
        ----------
        config : instance of ConfigParser
            the configuration, containing a [plot] section with options that
            control plotting

        nino34Index : xarray.dataArray
            nino34 timeseries to plot

        nino34Obs : xarray.dataArray
            nino34 observation

        nino3430 : xarray.dataArray
            subset of nino34 observations

        title : str
            the title of the plot

        obsTitle : str
            title of observational plot

        modelTitle : str
            title of model plot

        xLabel, yLabel : str
            axis labels

        fileout : str
            the file name to be written

        lineWidths : list of str
            control line width

        titleFontSize : int, optional
            the size of the title font

        figsize : tuple of floa  # {{{t, optional
            the size of the figure in inches

        dpi : int, optional
            the number of dots per inch of the figure, taken from section
            ``plot`` option ``dpi`` in the config file by default

        maxXTicks : int, optional
            the maximum number of tick marks that will be allowed along the x
            axis. This may need to be adjusted depending on the figure size and
            aspect ratio.

        Author
        ------
        Luke Van Roekel
        """
        if dpi is None:
            dpi = config.getint('plot', 'dpi')
        fig = plt.figure(figsize=figsize, dpi=dpi)

        if titleFontSize is None:
            titleFontSize = config.get('plot', 'titleFontSize')

        axis_font = {'size': config.get('plot', 'axisFontSize')}
        title_font = {'size': titleFontSize,
                      'color': config.get('plot', 'titleFontColor'),
                      'weight': config.get('plot', 'titleFontWeight')}
        if title is not None:
            fig.suptitle(title, y=0.92, **title_font)

        # Plot Nino34 Observation Time series
        plt.subplot(3, 1, 1)
        self._plot_nino_timeseries(plt, nino34Obs[2:-3].values,
                                   nino34Obs.Time[2:-3].values,
                                   xlabel, ylabel, obsTitle+' (Full Record)',
                                   calendar, axis_font, linewidths, maxXTicks)

        # Plot subset of the observational data set
        plt.subplot(3, 1, 2)
        self._plot_nino_timeseries(plt, nino3430.values, nino3430.Time.values,
                                   xlabel, ylabel, obsTitle+' (1976 - 2016)',
                                   calendar, axis_font, linewidths, maxXTicks)

        # Plot Nino34 model time series
        plt.subplot(3, 1, 3)
        self._plot_nino_timeseries(plt, nino34Index[2:-3].values,
                                   nino34Index.Time[2:-3].values,
                                   xlabel, ylabel, modelTitle, calendar,
                                   axis_font, linewidths, maxXTicks)
        minDays = nino34Index.Time[2:-3].values.min()
        maxDays = nino34Index.Time[2:-3].values.max()

        plot_xtick_format(plt, calendar, minDays, maxDays, maxXTicks)

        if fileout is not None:
            plt.savefig(fileout, dpi=dpi, bbox_inches='tight', pad_inches=0.1)

        if not config.getboolean('plot', 'displayToScreen'):
            plt.close()
        # }}}

    def _plot_nino_timeseries(self, plt, ninoIndex, time, xlabel, ylabel,
                              panelTitle, calendar, axis_font, linewidths,
                              maxXTicks):  # {{{
        '''
        Plot the nino time series on a subplot

        Parameters
        ----------
        ninoIndex : numpy.array
          nino34 Index values (can be obs or model)

        time : numpy.array
          time values for the nino index

        calendar : specified calendar for the plot

        maxXTicks : int, optional
            the maximum number of tick marks that will be allowed along the
            x axis. This may need to be adjusted depending on the figure size
            and aspect ratio.

        panelTitle : string
            string to label the subplot with

        xlabel : string
            string for x-axis label

        ylabel : string
            string for y-axis label

        Author
        ------
        Luke Van Roekel
        '''
        plt.title(panelTitle, y=1.06, **axis_font)
        y1 = ninoIndex
        nt = np.size(ninoIndex)

        y2 = np.zeros(nt)

        plt.plot(time, 0.4*np.ones(nt), '--k',
                 linewidth=linewidths)
        plt.plot(time, -0.4*np.ones(nt), '--k',
                 linewidth=linewidths)
        plt.fill_between(time, y1, y2, where=y1 > y2,
                         facecolor='red', interpolate=True, linewidth=0)
        plt.fill_between(time, y1, y2, where=y1 < y2,
                         facecolor='blue', interpolate=True, linewidth=0)

        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)
        # }}}

    def _write_xml(self, filePrefix, plotType):  # {{{
        caption = u'{} of El Niño 3.4 Climate Index'.format(plotType)
        write_image_xml(
            config=self.config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=u'El Niño 3.4 Climate Index',
            groupLink='nino34',
            thumbnailDescription=plotType,
            imageDescription=caption,
            imageCaption=caption)  # }}}

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
