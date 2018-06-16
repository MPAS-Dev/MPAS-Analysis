# -*- coding: utf-8 -*-
#

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import datetime
import xarray as xr
import pandas as pd
import numpy as np
from scipy import signal, stats
import matplotlib.pyplot as plt

from mpas_analysis.shared.climatology import climatology
from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.io.utility import build_config_full_path

from mpas_analysis.shared.timekeeping.utility import datetime_to_days, \
    string_to_days_since_date

from mpas_analysis.shared.io import open_mpas_dataset

from mpas_analysis.shared.plot.plotting import plot_xtick_format

from mpas_analysis.shared import AnalysisTask
from mpas_analysis.shared.html import write_image_xml


class IndexNino34(AnalysisTask):  # {{{
    '''
    A task for computing and plotting time series and spectra of the El Nino
    3.4 climate index

    Attributes
    ----------

    mpasTimeSeriesTask : ``MpasTimeSeriesTask``
        The task that extracts the time series from MPAS monthly output

    refConfig :  ``MpasAnalysisConfigParser``
        Configuration options for a reference run (if any)
    '''
    # Authors
    # -------
    # Luke Van Roekel, Xylar Asay-Davis

    def __init__(self, config, mpasTimeSeriesTask, refConfig=None):
        # {{{
        '''
        Construct the analysis task.

        Parameters
        ----------
        config :  ``MpasAnalysisConfigParser``
            Configuration options

        mpasTimeSeriesTask : ``MpasTimeSeriesTask``
            The task that extracts the time series from MPAS monthly output

        refConfig :  ``MpasAnalysisConfigParser``, optional
            Configuration options for a reference run (if any)
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call the constructor from the base class (AnalysisTask)
        super(IndexNino34, self).__init__(
            config=config,
            taskName='indexNino34',
            componentName='ocean',
            tags=['timeSeries', 'index', 'nino', 'publicObs'])

        self.mpasTimeSeriesTask = mpasTimeSeriesTask
        self.refConfig = refConfig

        self.run_after(mpasTimeSeriesTask)

        # }}}

    def setup_and_check(self):  # {{{
        '''
        Perform steps to set up the analysis and check for errors in the setup.
        '''
        # Authors
        # -------
        # Xylar Asay-Davis

        # first, call setup_and_check from the base class (AnalysisTask),
        # which will perform some common setup, including storing:
        #     self.runDirectory , self.historyDirectory, self.plotsDirectory,
        #     self.namelist, self.runStreams, self.historyStreams,
        #     self.calendar
        super(IndexNino34, self).setup_and_check()

        self.variableList = \
            ['timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']
        self.mpasTimeSeriesTask.add_variables(variableList=self.variableList)

        self.inputFile = self.mpasTimeSeriesTask.outputFile

        mainRunName = self.config.get('runs', 'mainRunName')

        config = self.config
        regionToPlot = config.get('indexNino34', 'region')

        if regionToPlot not in ['nino3.4', 'nino3', 'nino4']:
            raise ValueError('Unexpectes El Nino Index region {}'.format(
                    regionToPlot))
        ninoIndexNumber = regionToPlot[4:]

        self.xmlFileNames = []
        for filePrefix in ['nino{}_{}'.format(ninoIndexNumber, mainRunName),
                           'nino{}_spectra_{}'.format(ninoIndexNumber,
                                                      mainRunName)]:
            self.xmlFileNames.append('{}/{}.xml'.format(self.plotsDirectory,
                                                        filePrefix))

        # }}}

    def run_task(self):  # {{{
        '''
        Computes NINO34 index and plots the time series and power spectrum with
        95 and 99% confidence bounds
        '''
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        config = self.config
        calendar = self.calendar

        regionToPlot = config.get('indexNino34', 'region')

        ninoIndexNumber = regionToPlot[4:]

        self.logger.info("\nPlotting El Nino {} Index time series and power "
                         "spectrum....".format(ninoIndexNumber))

        self.logger.info('  Load SST data...')
        fieldName = 'nino'

        startDate = self.config.get('index', 'startDate')
        endDate = self.config.get('index', 'endDate')

        startYear = self.config.getint('index', 'startYear')
        endYear = self.config.getint('index', 'endYear')

        dataSource = config.get('indexNino34', 'observationData')

        observationsDirectory = build_config_full_path(
            config, 'oceanObservations', '{}Subdirectory'.format(fieldName))

        # specify obsTitle based on data path
        # These are the only data sets supported
        if dataSource == 'HADIsst':
            dataPath = "{}/HADIsst_nino34.nc".format(observationsDirectory)
            obsTitle = 'HADSST'
            refDate = '1870-01-01'
        elif dataSource == 'ERS_SSTv4':
            dataPath = "{}/ERS_SSTv4_nino34.nc".format(observationsDirectory)
            obsTitle = 'ERS SSTv4'
            refDate = '1800-01-01'
        else:
            raise ValueError('Bad value for config option observationData {} '
                             'in [indexNino34] section.'.format(dataSource))

        mainRunName = config.get('runs', 'mainRunName')

        # regionIndex should correspond to NINO34 in surface weighted Average
        # AM
        regions = config.getExpression('regions', 'regions')
        regionToPlot = config.get('indexNino34', 'region')
        regionIndex = regions.index(regionToPlot)

        # Load data:
        ds = open_mpas_dataset(fileName=self.inputFile,
                               calendar=calendar,
                               variableList=self.variableList,
                               startDate=startDate,
                               endDate=endDate)

        # Observations have been processed to the nino34Index prior to reading
        dsObs = xr.open_dataset(dataPath, decode_cf=False, decode_times=False)
        # add the days between 0001-01-01 and the refDate so we have a new
        # reference date of 0001-01-01 (like for the model Time)
        dsObs["Time"] = dsObs.Time + \
            string_to_days_since_date(dateString=refDate, calendar=calendar)
        nino34Obs = dsObs.sst

        self.logger.info('  Compute El Nino {} Index...'.format(
                ninoIndexNumber))
        varName = self.variableList[0]
        regionSST = ds[varName].isel(nOceanRegions=regionIndex)
        nino34Main = self._compute_nino34_index(regionSST, calendar)

        # Compute the observational index over the entire time range
        # nino34Obs = compute_nino34_index(dsObs.sst, calendar)

        self.logger.info(' Computing El Nino {} power spectra...'.format(
                ninoIndexNumber))
        spectraMain = self._compute_nino34_spectra(nino34Main)

        # Compute the observational spectra over the whole record
        spectraObs = self._compute_nino34_spectra(nino34Obs)

        # Compute the observational spectra over the last 30 years for
        # comparison. Only saving the spectra
        subsetEndYear = 2016
        if self.refConfig is None:
            subsetStartYear = 1976
        else:
            # make the subset the same length as the input data set
            subsetStartYear = subsetEndYear - (endYear - startYear)
        time_start = datetime_to_days(datetime.datetime(subsetStartYear, 1, 1),
                                      calendar=calendar)
        time_end = datetime_to_days(datetime.datetime(subsetEndYear, 12, 31),
                                    calendar=calendar)
        nino34Subset = nino34Obs.sel(Time=slice(time_start, time_end))
        spectraSubset = self._compute_nino34_spectra(nino34Subset)

        if self.refConfig is None:
            nino34s = [nino34Obs[2:-3], nino34Subset, nino34Main[2:-3]]
            titles = ['{} (Full Record)'.format(obsTitle),
                      '{} ({} - {})'.format(obsTitle, subsetStartYear,
                                            subsetEndYear),
                      mainRunName]
            spectra = [spectraObs, spectraSubset, spectraMain]
        else:
            baseDirectory = build_config_full_path(
                self.refConfig, 'output', 'timeSeriesSubdirectory')

            refFileName = '{}/{}.nc'.format(
                    baseDirectory, self.mpasTimeSeriesTask.fullTaskName)

            dsRef = open_mpas_dataset(
                    fileName=refFileName,
                    calendar=calendar,
                    variableList=self.variableList)

            regionSSTRef = dsRef[varName].isel(nOceanRegions=regionIndex)
            nino34Ref = self._compute_nino34_index(regionSSTRef, calendar)

            nino34s = [nino34Subset, nino34Main[2:-3], nino34Ref[2:-3]]
            refRunName = self.refConfig.get('runs', 'mainRunName')

            spectraRef = self._compute_nino34_spectra(nino34Ref)

            titles = ['{} ({} - {})'.format(obsTitle, subsetStartYear,
                                            subsetEndYear),
                      mainRunName,
                      'Ref: {}'.format(refRunName)]
            spectra = [spectraSubset, spectraMain, spectraRef]

        # Convert frequencies to period in years
        for s in spectra:
            s['period'] = \
                1.0 / (constants.eps + s['f']*constants.sec_per_year)

        self.logger.info(' Plot El Nino {} index and spectra...'.format(
                ninoIndexNumber))

        outFileName = '{}/nino{}_{}.png'.format(self.plotsDirectory,
                                                ninoIndexNumber, mainRunName)
        self._nino34_timeseries_plot(
                nino34s=nino34s,
                title=u'El Niño {} Index'.format(ninoIndexNumber),
                panelTitles=titles,
                outFileName=outFileName)

        self._write_xml(filePrefix='nino{}_{}'.format(ninoIndexNumber,
                                                      mainRunName),
                        plotType='Time Series',
                        ninoIndexNumber=ninoIndexNumber)

        outFileName = '{}/nino{}_spectra_{}.png'.format(self.plotsDirectory,
                                                        ninoIndexNumber,
                                                        mainRunName)
        self._nino34_spectra_plot(
                spectra=spectra,
                title=u'El Niño {} power spectrum'.format(ninoIndexNumber),
                panelTitles=titles,
                outFileName=outFileName)

        self._write_xml(filePrefix='nino{}_spectra_{}'.format(ninoIndexNumber,
                                                              mainRunName),
                        plotType='Spectra',
                        ninoIndexNumber=ninoIndexNumber)

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
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

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
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

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
        spectra = {'f': f, 'spectrum': pxxSmooth, 'conf99': mkov*scale*xHigh,
                   'conf95': mkov*scale*xLow, 'redNoise': mkov*scale}
        return spectra
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
        """
        # Authors
        # -------
        # Luke Van Roekel

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
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        nt = len(inputData)
        sp = (len(wgts) - 1) // 2
        runningMean = inputData.copy()
        for k in range(sp, nt-(sp+1)):
            runningMean[k] = sum(wgts*inputData[k-sp:k+sp+1].values)

        return runningMean  # }}}

    def _nino34_spectra_plot(self, spectra, title, panelTitles,
                             outFileName, lineWidth=2, xlabel='Period (years)',
                             ylabel=r'Power ($^o$C / cycles mo$^{-1}$)',
                             titleFontSize=None, figsize=(9, 21), dpi=None,
                             periodMin=1., periodMax=10.):
        # {{{
        """
        Plots the nino34 time series and power spectra in an image file
        Parameters
        ----------
        spectra : list of dict
            a dictionary for each panel returned from
            ``self._compute_nino34_spectra`` including entries
            ``period`` (periods to plot on x-axis), ``spectrum`` (nino34 power
            spectra), ``conf95`` (95% confidence level based on chi squared
            test), ``conf99`` (99% confidence level based on chi squared test)
            and ``redNoise`` (red noise fit to ``spectrum``)

        title : str
            the title of the plot

        panelTitles : list of str
            title of each panel of the plot

        outFileName : str
            the file name to be written

        lineWidth : int, optional
            control line width

        xLabel, yLabel : str, optional
            axis labels

        titleFontSize : int, optional
            the size of the title font

        figsize : tuple of float, optional
            the size of the figure in inches

        dpi : int, optional
            the number of dots per inch of the figure, taken from section
            ``plot`` option ``dpi`` in the config file by default

        periodMin, periodMax : float, optional
            the maximum and minimum periods (in years) to be plotted
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        config = self.config

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

        spectrumNames = ['spectrum', 'redNoise', 'conf95', 'conf99']
        colors = ['k', 'r', 'g', 'b']
        legends = ['Nino34 spectrum', 'Red noise fit',
                   '95% confidence threshold', '99% confidence threshold']

        maxYval = -1e20
        for plotIndex in range(3):
            x = spectra[plotIndex]['period']
            ys = [spectra[plotIndex][spectrumNames[curveIndex]] for curveIndex
                  in range(4)]
            maxYval = max(maxYval,
                          self._plot_size_y_axis(x=x, ys=ys, xmin=periodMin,
                                                 xmax=periodMax))

        for plotIndex in range(3):
            plt.subplot(3, 1, plotIndex+1)

            period = spectra[plotIndex]['period']
            for curveIndex in range(4):
                spectrum = spectra[plotIndex][spectrumNames[curveIndex]]
                plt.plot(period[2:-3], spectrum[2:-3], colors[curveIndex],
                         linewidth=lineWidth, label=legends[curveIndex])
            plt.xlim(10, 1)

            plt.legend(loc='upper right')
            plt.ylim(0, 0.9*maxYval)

            if panelTitles[plotIndex] is not None:
                plt.title(panelTitles[plotIndex], **title_font)
            if xlabel is not None:
                plt.xlabel(xlabel, **axis_font)
            if ylabel is not None:
                plt.ylabel(ylabel, **axis_font)

        plt.tight_layout(rect=[0, 0.03, 1, 0.90])

        if outFileName is not None:
            fig.savefig(outFileName, dpi=dpi, bbox_inches='tight',
                        pad_inches=0.1)

        plt.close()
        # }}}

    def _nino34_timeseries_plot(self, nino34s, title, panelTitles, outFileName,
                                xlabel='Time (years)', ylabel='($\degree$C)',
                                titleFontSize=None, figsize=(9, 21), dpi=None,
                                maxXTicks=20, lineWidth=2):
        # {{{
        """
        Plots the nino34 time series and power spectra in an image file

        Parameters
        ----------
        nino34s : list of xarray.dataArray
            nino34 timeseries to plot in each panel

        title : str
            the title of the plot

        panelTitles : list of str
            title of each panel of the plot

        outFileName : str
            the file name to be written

        xLabel, yLabel : str
            axis labels

        titleFontSize : int, optional
            the size of the title font

        figsize : tuple of float, optional
            the size of the figure in inches

        dpi : int, optional
            the number of dots per inch of the figure, taken from section
            ``plot`` option ``dpi`` in the config file by default

        lineWidth : int, optional
            control line width
        """
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        config = self.config
        calendar = self.calendar

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

        for plotIndex in range(3):
            plt.subplot(3, 1, plotIndex+1)
            index = nino34s[plotIndex].values
            time = nino34s[plotIndex].Time.values
            self._plot_nino_timeseries(index, time, xlabel, ylabel,
                                       panelTitles[plotIndex],
                                       title_font, axis_font, lineWidth)

            minDays = time.min()
            maxDays = time.max()

            plot_xtick_format(calendar, minDays, maxDays, maxXTicks)

        plt.tight_layout(rect=[0, 0.03, 1, 0.90])

        if outFileName is not None:
            plt.savefig(outFileName, dpi=dpi, bbox_inches='tight',
                        pad_inches=0.1)

        plt.close()
        # }}}

    def _plot_nino_timeseries(self, ninoIndex, time, xlabel, ylabel,
                              panelTitle, title_font, axis_font,
                              lineWidth):  # {{{
        '''
        Plot the nino time series on a subplot

        Parameters
        ----------
        ninoIndex : numpy.array
          nino34 Index values (can be obs or model)

        time : numpy.array
          time values for the nino index

        xlabel : string
            string for x-axis label

        ylabel : string
            string for y-axis label

        panelTitle : string
            string to label the subplot with

        lineWidth : list of str
            control line width
        '''
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        plt.title(panelTitle, y=1.06, **title_font)
        y1 = ninoIndex
        nt = np.size(ninoIndex)

        y2 = np.zeros(nt)

        plt.plot(time, 0.4*np.ones(nt), '--k',
                 linewidth=lineWidth)
        plt.plot(time, -0.4*np.ones(nt), '--k',
                 linewidth=lineWidth)
        plt.fill_between(time, y1, y2, where=y1 > y2,
                         facecolor='red', interpolate=True, linewidth=0)
        plt.fill_between(time, y1, y2, where=y1 < y2,
                         facecolor='blue', interpolate=True, linewidth=0)

        if xlabel is not None:
            plt.xlabel(xlabel, **axis_font)
        if ylabel is not None:
            plt.ylabel(ylabel, **axis_font)
        # }}}

    def _write_xml(self, filePrefix, plotType, ninoIndexNumber):  # {{{
        caption = u'{} of El Niño {} Climate Index'.format(plotType,
                                                           ninoIndexNumber)
        write_image_xml(
            config=self.config,
            filePrefix=filePrefix,
            componentName='Ocean',
            componentSubdirectory='ocean',
            galleryGroup=u'El Niño {} Climate Index'.format(ninoIndexNumber),
            groupLink='nino',
            thumbnailDescription=plotType,
            imageDescription=caption,
            imageCaption=caption)  # }}}

    def _plot_size_y_axis(self, x, ys, xmin, xmax):
        '''
        Get the maximum y value over the given range of x values

        Parameters
        ----------
        x : numpy.array
           x values

        ys : list of numpy.array
           a list of curves (y values)

        xmin : float
            The minimum x value

        xmax : float, optional
            The maximum x values
        '''
        # Authors
        # -------
        # Luke Van Roekel, Xylar Asay-Davis

        mask = np.logical_and(x >= xmin, x <= xmax)

        # find maximum value of three curves plotted
        maxY = -1E20
        for y in ys:
            maxY = max(y[mask].max(), maxY)
            # check the function interpolated to the max/min as well
            # Note: flipping the axis so x is in increasing order
            maxY = max(np.interp(xmin, x[::-1], y[::-1]), maxY)
            maxY = max(np.interp(xmax, x[::-1], y[::-1]), maxY)

        return maxY

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
