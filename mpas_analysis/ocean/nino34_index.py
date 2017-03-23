"""
Computes NINO34 index and plots the time series and power spectra

Author
------
Luke Van Roekel

Last Modified
-------------
03/21/2017
"""
from ..shared.plot.plotting import nino34_timeseries_plot, nino34_spectra_plot
import xarray as xr
import pandas as pd
import numpy as np

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath
from ..shared.climatology import climatology
from ..shared.constants import constants

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time


def nino34_index(config, streamMap=None, variableMap=None):
    # {{{
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
    Luke Van Roekel

    Last Modified
    -------------
    03/22/2017
    """

    # Define/read in general variables
    print '  Load SST data...'
    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    streamsFileName = config.get('input', 'oceanStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    calendar = namelist.get('config_calendar_type')
    simulationStartTime = get_simulation_start_time(streams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('index', 'startDate')
    endDate = config.get('index', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    fileNames = streams.readpath(streamName, startDate=startDate,
                                 endDate=endDate,  calendar=calendar)
    print 'Reading files {} through {}'.format(fileNames[0], fileNames[-1])

    mainRunName = config.get('runs', 'mainRunName')
    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')

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

    SSTregions = ds.avgSurfaceTemperature
    print '  Compute NINO3.4 index...'
    nino34 = compute_nino34_index(SSTregions[:, regionIndex], config)

    print ' Computing NINO3.4 power spectra...'
    f, spectra, conf99, conf95, redNoise = compute_nino34_spectra(nino34,
                                                                  config)

    # Convert frequencies to period in years
    f = 1.0 / (constants.eps + f*constants.sec_per_year)
    print ' Plot NINO3.4 index and spectra...'

    figureName = '{}/NINO34_{}.png'.format(plotsDirectory, mainRunName)
    nino34_timeseries_plot(config, nino34, 'NINO 3.4 Index',
                           figureName, linewidths=2, calendar=calendar)

    figureName = '{}/NINO34_spectra_{}.png'.format(plotsDirectory, mainRunName)
    nino34_spectra_plot(config, f, spectra, conf95, conf99, redNoise,
                        'NINO3.4 power spectrum', figureName, linewidths=2)


# }}}
def compute_nino34_index(regionSST, config):
    # {{{
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
    regionSST : xarray.dataArray object
       values of SST in the nino region

    config : MpasConfigParser object
        the config options

    Returns
    -------
    xarray.DataArray object containing the nino34index

    Author
    ------
    Luke Van Roekel

    Last Modified
    -------------
    03/22/2017
    """

    assert isinstance(regionSST, xr.core.dataarray.DataArray)

    inDirectory = config.get('input', 'baseDirectory')
    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)
    calendar = namelist.get('config_calendar_type')

    # Compute monthly average and anomaly of climatology of SST
    monthlyClimatology = climatology.compute_monthly_climatology(regionSST,
                                                                 calendar)
    anomalySST = regionSST.groupby('month') - monthlyClimatology

    return _running_mean(anomalySST.to_pandas())


# }}}
def compute_nino34_spectra(nino34Index, config):
    # {{{
    """
    Computes power spectra of Nino34 index.

    nino34Index is the NINO index computed by compute_nino34_index

    The algorithm follows the NCL cvdp package see
    http://www.cesm.ucar.edu/working_groups/CVC/cvdp/code.html

    Parameters
    ----------
    nino34Index : xarray.DataArray object
        nino34Index for analysis

    config : instance of the MPAS configParser

    Returns
    -------
    pxxSmooth : xarray.DataArray object
        nino34Index power spectra that has been smoothed with a 5-point
        running mean

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
    Luke Van Roekel

    Last Modified
    -------------
    03/22/2017
    """

    from scipy import signal, stats

    detrendedNino34 = signal.detrend(nino34Index.values)

    window = signal.tukey(len(detrendedNino34), alpha=0.1)
    f, Pxx = signal.periodogram(window * detrendedNino34,
                                1.0 / constants.sec_per_month,
                                scaling='spectrum')

    # computes power spectra, smoothed with 5 point running mean
    pxxSmooth = _running_mean(pd.Series(Pxx))

    # compute 99 and 95% confidence intervals and red-noise process
    # Uses Chi squared test

    r = _autocorr(detrendedNino34)[0, 1]
    r2 = 2.*r
    rsq = r**2

    # In the temp2 variable, f is converted to give wavenumber, i.e.
    # 0,1,2,...,N/2
    temp2 = r2*np.cos(2.*np.pi*f*constants.sec_per_month)
    mkov = 1. / (1. + rsq - temp2)

    sum1 = np.sum(mkov)
    sum2 = np.sum(pxxSmooth.values)
    scale = sum2 / sum1

    xLow = stats.chi2.interval(0.95, 4)[1]/4.
    xHigh = stats.chi2.interval(0.99, 4)[1]/4.

    # return Spectra, 99% confidence level, 95% confidence level,
    #        and Red-noise fit
    return f, pxxSmooth, mkov*scale*xHigh, mkov*scale*xLow, mkov*scale


# }}}
def _autocorr(x, t=1):
    # {{{
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

    return np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))


# }}}
def _running_mean(inputData):
    # {{{
    """
    Calculates 5-month running mean for NINO index and the spectra

    Author
    ------
    Luke Van Roekel

    Last Modified
    -------------
    03/22/2017
    """

    runningMean = pd.Series.rolling(inputData, 5, center=True, min_periods=1).mean()
    return xr.DataArray.from_series(runningMean)


# }}}
# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
