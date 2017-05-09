import datetime
import xarray as xr
import pandas as pd
import numpy as np
from scipy import signal, stats
import os

from ..shared.climatology import Climatology
from ..shared.constants import constants
from ..shared.io import build_config_full_path
from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    add_years_months_days_in_month

from ..shared.plot.plotting import nino34_timeseries_plot, nino34_spectra_plot

from ..shared.analysis_task import AnalysisTask


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
        #     self.calendar, self.namelistMap, self.streamMap, self.variableMap
        super(IndexNino34, self).setup_and_check()

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

        # get a list of timeSeriesStats output files from the streams file,
        # reading only those that are between the start and end dates
        startDate = config.get('index', 'startDate')
        endDate = config.get('index', 'endDate')
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

        streamName = self.historyStreams.find_stream(
            self.streamMap['timeSeriesStats'])
        fileNames = self.historyStreams.readpath(streamName,
                                                 startDate=startDate,
                                                 endDate=endDate,
                                                 calendar=calendar)
        print '\n  Reading files:\n' \
              '    {} through\n    {}'.format(
                  os.path.basename(fileNames[0]),
                  os.path.basename(fileNames[-1]))
        mainRunName = config.get('runs', 'mainRunName')

        # regionIndex should correspond to NINO34 in surface weighted Average
        # AM
        regionIndex = config.getint('indexNino34', 'regionIndicesToPlot')

        # Load data:
        varList = ['avgSurfaceTemperature']
        ds = open_multifile_dataset(fileNames=fileNames,
                                    calendar=calendar,
                                    config=config,
                                    simulationStartTime=simulationStartTime,
                                    timeVariableName='Time',
                                    variableList=varList,
                                    variableMap=self.variableMap,
                                    startDate=startDate,
                                    endDate=endDate)

        # Observations have been processed to the nino34Index prior to reading
        dsObs = xr.open_dataset(dataPath)
        nino34Obs = dsObs.sst

        print '  Compute NINO3.4 index...'
        regionSST = ds.avgSurfaceTemperature.isel(nOceanRegions=regionIndex)
        nino34 = self._compute_nino34_index(regionSST)

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
        nino34_timeseries_plot(config, nino34, nino34Obs, nino3430,
                               'NINO 3.4 Index', modelTitle, obsTitle,
                               figureName, linewidths=2, calendar=calendar)

        figureName = '{}/NINO34_spectra_{}.png'.format(self.plotsDirectory,
                                                       mainRunName)
        nino34_spectra_plot(config, f, spectra, conf95, conf99, redNoise,
                            fObs, f30, spectraObs, conf95Obs, conf99Obs,
                            redNoiseObs, spectra30yrs, conf9530, conf9930,
                            redNoise30, 'NINO3.4 power spectrum', modelTitle,
                            obsTitle, figureName, linewidths=2)
    # }}}

    def _compute_nino34_index(self, regionSST):  # {{{
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
        regionSST = add_years_months_days_in_month(regionSST, self.calendar)

        # Compute monthly average and anomaly of climatology of SST
        monthlyClimatology = Climatology(task=self)
        monthlyClimatology.compute_monthly(regionSST, maskVaries=False)

        anomaly = regionSST.groupby('month') - monthlyClimatology.dataSet

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

# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
