from ..shared.plot.plotting import nino34_timeseries_plot, nino34_spectra_plot
import pandas as pd
import numpy as np
from scipy import signal, stats

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath
from ..shared.climatology import climatology

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time, \
    date_to_days, days_to_datetime


def nino34_timeseries(config, streamMap=None, variableMap=None, SSTregions=None):
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
    03/21/2017
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
    startDate = config.get('timeSeries', 'startDate')
    endDate = config.get('timeSeries', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    fileNames = streams.readpath(streamName, startDate=startDate,
                                 endDate=endDate,  calendar=calendar)
    print 'Reading files {} through {}'.format(fileNames[0], fileNames[-1])

    mainRunName = config.get('runs', 'mainRunName')
    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')

    plotTitles = config.getExpression('regions', 'plotTitles')
#   regionIndex should correspond to NINO34 in surface weighted Average AM
    regionIndex = config.getint('timeSeriesNino34', 'regionIndex')

#   Load data:

    varList = ['avgSurfaceTemperature']
    ds = open_multifile_dataset(fileNames=fileNames,
                                calendar=calendar,
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
#   NOTE: first 2 and last two values ommitted due to rolling mean
    f, spectra, conf99, conf95, redNoise = compute_nino34_spectra(nino34, config)

    start = days_to_datetime(np.amin(ds.Time.min()), calendar=calendar)
    end = days_to_datetime(np.amax(ds.Time.max()), calendar=calendar)

#   Convert frequencies to period in years
    f = 1.0 / (1E-9 + f*86400.*365)
    print ' Plot NINO3.4 index and spectra...'
    xLabel = 'Time [years]'
    yLabel = '[$^\circ$ C]'

    figureName = '{}/NINO34_{}.png'.format(plotsDirectory, mainRunName)
    nino34_timeseries_plot(config, nino34, 'NINO 3.4 Index',
    			   xLabel, yLabel, figureName, linewidths=3, calendar=calendar)

    xLabel = 'Period (years)'
    yLabel = 'Power Spectrum'
    figureName = '{}/NINO34_spectra_{}.png'.format(plotsDirectory, mainRunName)
    nino34_spectra_plot(config, f, spectra, conf95, conf99, redNoise,
			'NINO3.4 power spectrum', xLabel,
			yLabel, figureName, linewidths=3)

 
def compute_nino34_index(regionSST, config):
    """
    Computes nino34 index time series.  It follow the standard nino34 algorithm, i.e.,

    1) Compute monthly average SST in the region
    2) Computes anomalous SST
    3) Performs a 5 point running mean over the anomalies

    This routine requires regionSST to be the SSTs in the nino3.4 region ONLY.
    It is defined as lat > -5S and lat < 5N and lon > 190E and lon < 240E.

    Further, regionSST must be a dataArray

    Author: Luke Van Roekel
    Last Modified: 03/20/2017
    """

    inDirectory = config.get('input', 'baseDirectory')
    namelistFileName = config.get('input', 'oceanNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)
    calendar = namelist.get('config_calendar_type')

    # Compute monthly average climatology of SST
    monthlyClimatology = climatology.compute_monthly_climatology(regionSST, calendar)

    anomalySST = regionSST.groupby('month') - monthlyClimatology

    return _running_mean(anomalySST.to_pandas())


def compute_nino34_spectra(nino34Index, config):
    """
    Computes power spectra of Nino34 index.

    nino34Index is the NINO index computed by compute_nino34_index

    The algorithm follows the NCL cvdp package see
    http://www.cesm.ucar.edu/working_groups/CVC/cvdp/code.html

    NOTE: this routine requires the signal and stats module from scipy

    Author: Luke Van Roekel
    Last Modified: 03/21/2017
    """

    detrendedNino34 = signal.detrend(nino34Index.values)

    f, Pxx = signal.periodogram(nino34Index.values, 1.0 / (86400.*30.), 
				window='hamming', scaling='spectrum')

#   computes power spectra, smoothed with 5 point running mean
    pxxSmooth = _running_mean(pd.Series(Pxx))
    
#   compute 99 and 95% confidence intervals and red-noise process
#   Uses Chi squared test

    r = _autocorr(detrendedNino34)[0, 1]
    r2 = 2.*r
    rsq = r**2

#   In the temp2 variable, f is converted to give wavenumber i.e., 0,1,2,...,N/2
    temp2 = r2*np.cos(2.*np.pi*f*86400.*30.)
    mkov = 1. / (1. + rsq - temp2)

    sum1 = np.sum(mkov)
    sum2 = np.sum(pxxSmooth)
    scale = sum2 / sum1

    xLow = stats.chi2.interval(0.95, 2)[1]/2.
    xHigh = stats.chi2.interval(0.99, 2)[1]/2.

#   return Spectra, 99% confidence level, 95% confidence level, and Red-noise fit
    return f, pxxSmooth, mkov*scale*xHigh, mkov*scale*xLow, mkov*scale


def _autocorr(x, t=1):
    """
    Computes lag one auto-correlation for the NINO34 spectra calculation

    Author: Luke Van Roekel
    Last Modified: 03/20/2017
    """

    return np.corrcoef(np.array([x[0:len(x)-t], x[t:len(x)]]))


def _running_mean(inputData):
    """
    Calculates 5-point running mean for NINO index and the spectra

    Author: Luke Van Roekel
    Last Modified: 03/20/2017
    """

    return pd.Series.rolling(inputData, 5, center=True, min_periods=1).mean()
