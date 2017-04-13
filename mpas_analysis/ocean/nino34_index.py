"""
Computes NINO34 index and plots the time series and power spectra

Author
------
Luke Van Roekel, Xylar Asay-Davis

Last Modified
-------------
04/10/2017
"""

import datetime
import xarray as xr
import pandas as pd
import numpy as np
from scipy import signal, stats
import os

from ..shared.climatology import climatology
from ..shared.constants import constants
from ..shared.io.utility import build_config_full_path
from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.plot.plotting import nino34_timeseries_plot, nino34_spectra_plot

from ..shared.analysis_task import setup_task


def nino34_index(config, streamMap=None, variableMap=None):  # {{{
    """
    Computes NINO34 index and plots the time series and power spectrum with
    95 and 99% confidence bounds

    Parameters
    ----------
    config: Instance of MpasAnalysisConfigParser containing configuration
    options.

    streamMap: dict, optional
        a dictionary of MPAS-O variable names that map to their
        mpas_analysis counterparts.

    variableMap : dict, optional
        a dictionary of MPAS-O variable names that map
        to their mpas_analysis counterparts.

    Author
    ------
    Luke Van Roekel, Xylar Asay-Davis

    Last Modified
    -------------
    04/10/2017
    """

    print '  Load SST data...'
    field = 'nino'

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    simulationStartTime = get_simulation_start_time(runStreams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('index', 'startDate')
    endDate = config.get('index', 'endDate')
    dataSource = config.get('indexNino34', 'observationData')

    observationsDirectory = build_config_full_path(
        config, 'oceanObservations', '{}Subdirectory'.format(field))

    # specify obsTitle based on data path
    # These are the only data sets supported
    if dataSource == 'HADIsst':
        dataPath = "{}/HADIsst_nino34.nc".format(observationsDirectory)
        obsTitle = 'HADSST'
    else:
        dataPath = "{}/ERS_SSTv4_nino34.nc".format(observationsDirectory)
        obsTitle = 'ERS SSTv4'

    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    fileNames = historyStreams.readpath(streamName, startDate=startDate,
                                        endDate=endDate,  calendar=calendar)
    print '\n  Reading files:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(fileNames[0]),
              os.path.basename(fileNames[-1]))
    mainRunName = config.get('runs', 'mainRunName')

    # regionIndex should correspond to NINO34 in surface weighted Average AM
    regionIndex = config.getint('indexNino34', 'regionIndicesToPlot')

    # Load data:
    varList = ['avgSurfaceTemperature']
    ds = open_multifile_dataset(fileNames=fileNames,
                                calendar=calendar,
                                config=config,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=varList,
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    # Observations have been processed to the nino34Index prior to reading
    dsObs = xr.open_dataset(dataPath)
    nino34Obs = dsObs.sst

    print '  Compute NINO3.4 index...'
    regionSST = ds.avgSurfaceTemperature.isel(nOceanRegions=regionIndex)
    nino34 = compute_nino34_index(regionSST, calendar)

    # Compute the observational index over the entire time range
#    nino34Obs = compute_nino34_index(dsObs.sst, calendar)

    print ' Computing NINO3.4 power spectra...'
    f, spectra, conf99, conf95, redNoise = compute_nino34_spectra(nino34)

    # Compute the observational spectra over the whole record
    fObs, spectraObs, conf99Obs, conf95Obs, redNoiseObs = compute_nino34_spectra(nino34Obs)

    # Compute the observational spectra over the last 30 years for comparison
    # Only saving the spectra
    time_start = datetime.datetime(1976, 1, 1)
    time_end = datetime.datetime(2016, 12, 31)
    nino3430 = nino34Obs.sel(Time=slice(time_start, time_end))
    f30, spectra30yrs, conf9930, conf9530, redNoise30 = compute_nino34_spectra(nino3430)

    # Convert frequencies to period in years
    f = 1.0 / (constants.eps + f*constants.sec_per_year)
    fObs = 1.0 / (constants.eps + fObs*constants.sec_per_year)
    f30 = 1.0 / (constants.eps + f30*constants.sec_per_year)

    print ' Plot NINO3.4 index and spectra...'

    figureName = '{}/NINO34_{}.png'.format(plotsDirectory, mainRunName)
    modelTitle = "{}".format(mainRunName)
    nino34_timeseries_plot(config, nino34, nino34Obs, nino3430, 'NINO 3.4 Index',
                           modelTitle, obsTitle, figureName, linewidths=2,
                           calendar=calendar)

    figureName = '{}/NINO34_spectra_{}.png'.format(plotsDirectory, mainRunName)
    nino34_spectra_plot(config, f, spectra, conf95, conf99, redNoise,
                        fObs, f30, spectraObs, conf95Obs, conf99Obs, redNoiseObs,
                        spectra30yrs, conf9530, conf9930, redNoise30,
                        'NINO3.4 power spectrum', modelTitle,
                        obsTitle, figureName, linewidths=2)
    # }}}


def compute_nino34_index(regionSST, calendar):  # {{{
    """
    Computes nino34 index time series.  It follow the standard nino34
    algorithm, i.e.,

      1. Compute monthly average SST in the region
      2. Computes anomalous SST
      3. Performs a 5 month running mean over the anomalies

    This routine requires regionSST to be the SSTs in the nino3.4 region ONLY.
    It is defined as lat > -5S and lat < 5N and lon > 190E and lon < 240E.

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

    Last Modified
    -------------
    04/08/2017
    """

    if not isinstance(regionSST, xr.core.dataarray.DataArray):
        raise ValueError('regionSST should be an xarray DataArray')

    # add 'month' data array so we can group by month below.
    regionSST = climatology.add_years_months_days_in_month(regionSST, calendar)

    # Compute monthly average and anomaly of climatology of SST
    monthlyClimatology = \
        climatology.compute_monthly_climatology(regionSST, maskVaries=False)

    anomaly = regionSST.groupby('month') - monthlyClimatology

    # Remove the long term trend from the anomalies
    detrendedAnomal = signal.detrend(anomaly.values)
    anomaly.values = detrendedAnomal

    # Compute 5 month running mean
    wgts = np.ones(5) / 5.
    return _running_mean(anomaly, wgts)  # }}}


def compute_nino34_spectra(nino34Index):  # {{{
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

    Last Modified
    -------------
    04/10/2017
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

    pxxSmooth = _running_mean(pd.Series(Pxx), wgts) / constants.sec_per_month

    # compute 99 and 95% confidence intervals and red-noise process
    # Uses Chi squared test

    r = _autocorr(ninoIndex)[0, 1]
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
    return f, pxxSmooth, mkov*scale*xHigh, mkov*scale*xLow, mkov*scale  # }}}


def _autocorr(x, t=1):  # {{{
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

    Last Modified
    -------------
    03/22/2017
    """

    return np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))  # }}}


def _running_mean(inputData, wgts):  # {{{
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

    Last Modified
    -------------
    04/10/2017
    """

    nt = len(inputData)
    sp = (len(wgts) - 1)/2
    runningMean = inputData.copy()
    for k in range(sp, nt-(sp+1)):
        runningMean[k] = sum(wgts*inputData[k-sp:k+sp+1].values)

    return runningMean  # }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
