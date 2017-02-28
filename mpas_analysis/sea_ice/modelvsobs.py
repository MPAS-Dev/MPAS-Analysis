"""
General comparison of 2-d model fields against data.  Currently only supports
sea ice concentration (sic) and sea ice thickness (sit)

Author: Xylar Asay-Davis, Milena Veneziani
Last Modified: 02/26/2017
"""

import os
import os.path
import matplotlib.pyplot as plt
import matplotlib.colors as cols

import numpy.ma as ma

import netCDF4

from ..shared.constants import constants

from ..shared.climatology import climatology

from ..shared.plot.plotting import plot_polar_comparison

from ..shared.io import NameList, StreamsFile
from ..shared.io.utility import buildConfigFullPath

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time


def seaice_modelvsobs(config, streamMap=None, variableMap=None):
    """
    Performs analysis of sea-ice properties by comparing with
    previous model results and/or observations.

    config is an instance of MpasAnalysisConfigParser containing configuration
    options.

    If present, streamMap is a dictionary of MPAS-O stream names that map to
    their mpas_analysis counterparts.

    If present, variableMap is a dictionary of MPAS-O variable names that map
    to their mpas_analysis counterparts.

    Author: Xylar Asay-Davis, Milena Veneziani
    Last Modified: 02/26/2017
    """

    # read parameters from config file
    inDirectory = config.get('input', 'baseDirectory')

    namelistFileName = config.get('input', 'seaIceNamelistFileName')
    namelist = NameList(namelistFileName, path=inDirectory)

    streamsFileName = config.get('input', 'seaIceStreamsFileName')
    streams = StreamsFile(streamsFileName, streamsdir=inDirectory)

    calendar = namelist.get('config_calendar_type')
    try:
        simulationStartTime = get_simulation_start_time(streams)
    except IOError:
        # try the ocean stream instead
        oceanStreamsFileName = config.get('input', 'oceanStreamsFileName')
        oceanStreams = StreamsFile(oceanStreamsFileName, streamsdir=inDirectory)
        simulationStartTime = get_simulation_start_time(oceanStreams)
        oceanStreams.close()
    # get a list of timeSeriesStatsMonthly output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = streams.find_stream(streamMap['timeSeriesStats'])
    fileNames = streams.readpath(streamName, startDate=startDate,
                                 endDate=endDate, calendar=calendar)
    print 'Reading files {} through {}'.format(fileNames[0], fileNames[-1])

    try:
        restartFileName = streams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-Sea Ice restart file found: need at least one '
                      'restart file for modelvsobs calculation')

    # Load data
    print "  Load sea-ice data..."
    ds = open_multifile_dataset(fileNames=fileNames,
                                calendar=calendar,
                                simulationStartTime=simulationStartTime,
                                timeVariableName='Time',
                                variableList=['iceAreaCell',
                                              'iceVolumeCell'],
                                variableMap=variableMap,
                                startDate=startDate,
                                endDate=endDate)

    # Compute climatologies (first motnhly and then seasonally)
    print "  Compute seasonal climatologies..."

    monthlyClimatology = climatology.compute_monthly_climatology(ds, calendar)

    mpasInterolationData = \
        climatology.init_model_interpolation(restartFileName)

    _compute_and_plot_concentration(config, monthlyClimatology,
                                    mpasInterolationData)

    _compute_and_plot_thickness(config, monthlyClimatology,
                                mpasInterolationData)


def _compute_and_plot_concentration(config, monthlyClimatology,
                                    mpasInterolationData):
    """
    Given a config file, monthly climatology on the mpas grid, and the data
    necessary to perform horizontal interpolation to a comparison grid,
    computes seasonal climatologies and plots model results, observations
    and biases in sea-ice concentration.

    Parameters
    ----------
    config : an instance of MpasConfigParser

    monthlyClimatology : an xarray data set containing a monthly climatology

    mpasInterolationData : a tuple returned by init_model_interpolation that
        contains the information needed to interpolate from the MPAS mesh
        to the comparison grid

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    02/26/2017
    """

    print "  Make ice concentration plots..."

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')
    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')

    subtitle = "Ice concentration"

    hemisphereSeasons = {'JFM': ('NH', 'Winter'),
                         'JAS': ('NH', 'Summer'),
                         'DJF': ('SH', 'Winter'),
                         'JJA': ('SH', 'Summer')}

    lonTarg = mpasInterolationData[2]
    latTarg = mpasInterolationData[3]

    for months in hemisphereSeasons:
        hemisphere, season = hemisphereSeasons[months]
        monthsValue = constants.monthDictionary[months]

        iceConcentration = climatology.interpolate_model_climatology(
            mpasInterolationData, monthlyClimatology, "iceAreaCell",
            monthsValue)

        if hemisphere == 'NH':
            plotProjection = 'npstere'
        else:
            plotProjection = 'spstere'

        resultConcContourValues = config.getExpression(
            'regriddedSeaIceConcThick',
            'resultConc{}ContourValues'.format(season))
        resultColormap = plt.get_cmap(config.get(
            'regriddedSeaIceConcThick', 'resultColormap'))
        resultColormapIndices = config.getExpression(
            'regriddedSeaIceConcThick', 'resultColormapIndices')
        resultColormap = cols.ListedColormap(
            resultColormap(resultColormapIndices), "resultColormap")

        differenceConcContourValues = config.getExpression(
            'regriddedSeaIceConcThick',
            'differenceConc{}ContourValues'.format(season))
        differenceColormap = plt.get_cmap(config.get(
            'regriddedSeaIceConcThick', 'differenceColormap'))
        differenceColormapIndices = config.getExpression(
            'regriddedSeaIceConcThick', 'differenceColormapIndices')
        differenceColormap = cols.ListedColormap(
            differenceColormap(differenceColormapIndices),
            "differenceColormap")

        referenceLongitude = config.getfloat(
            'regriddedSeaIceConcThick',
            'referenceLongitude{}'.format(hemisphere))
        minimumLatitude = config.getfloat(
            'regriddedSeaIceConcThick',
            'minimumLatitude{}'.format(hemisphere))

        # ice concentrations from NASATeam (or Bootstrap) algorithm
        for obsName in ['NASATeam', 'Bootstrap']:

            obsFileNames = buildConfigFullPath(
                config, 'seaIceObservations',
                'concentration{}{}_{}'.format(obsName, hemisphere, months))

            if not os.path.isfile(obsFileNames):
                raise SystemExit("Obs file {} not found. Exiting...".format(
                    obsFileNames))

            ncFile = netCDF4.Dataset(obsFileNames, mode='r')
            obsIceConcentration = ncFile.variables["AICE"][:].T
            ncFile.close()

            difference = iceConcentration - obsIceConcentration

            title = "{} ({}, years {:04d}-{:04d})".format(
                subtitle, months, startYear, endYear)
            fileout = "{}/iceconc{}{}_{}_{}_years{:04d}-{:04d}.png".format(
                plotsDirectory, obsName, hemisphere, mainRunName,
                months, startYear, endYear)
            plot_polar_comparison(
                config,
                lonTarg,
                latTarg,
                iceConcentration,
                obsIceConcentration,
                difference,
                resultColormap,
                resultConcContourValues,
                differenceColormap,
                differenceConcContourValues,
                title=title,
                fileout=fileout,
                plotProjection=plotProjection,
                latmin=minimumLatitude,
                lon0=referenceLongitude,
                modelTitle=mainRunName,
                obsTitle="Observations (SSM/I {})".format(obsName),
                diffTitle="Model-Observations",
                cbarlabel="fraction")


def _compute_and_plot_thickness(config, monthlyClimatology,
                                mpasInterolationData):
    """
    Given a config file, monthly climatology on the mpas grid, and the data
    necessary to perform horizontal interpolation to a comparison grid,
    computes seasonal climatologies and plots model results, observations
    and biases in sea-ice thickness.

    Parameters
    ----------
    config : an instance of MpasConfigParser

    monthlyClimatology : an xarray data set containing a monthly climatology

    mpasInterolationData : a tuple returned by init_model_interpolation that
        contains the information needed to interpolate from the MPAS mesh
        to the comparison grid

    Authors
    -------
    Xylar Asay-Davis, Milena Veneziani

    Last Modified
    -------------
    02/26/2017
    """

    print "  Make ice thickness plots..."

    subtitle = "Ice thickness"

    plotsDirectory = buildConfigFullPath(config, 'output', 'plotsSubdirectory')
    mainRunName = config.get('runs', 'mainRunName')
    startYear = config.getint('climatology', 'startYear')
    endYear = config.getint('climatology', 'endYear')

    lonTarg = mpasInterolationData[2]
    latTarg = mpasInterolationData[3]

    for months in ['FM', 'ON']:
        monthsValue = constants.monthDictionary[months]
        iceThickness = climatology.interpolate_model_climatology(
            mpasInterolationData, monthlyClimatology, "iceVolumeCell",
            monthsValue)

        for hemisphere in ['NH', 'SH']:

            resultThickContourValues = config.getExpression(
                'regriddedSeaIceConcThick',
                'resultThick{}ContourValues'.format(hemisphere))
            resultColormap = plt.get_cmap(config.get(
                'regriddedSeaIceConcThick', 'resultColormap'))
            resultColormapIndices = config.getExpression(
                'regriddedSeaIceConcThick', 'resultColormapIndices')
            resultColormap = cols.ListedColormap(
                resultColormap(resultColormapIndices), "resultColormap")

            differenceThickContourValues = config.getExpression(
                'regriddedSeaIceConcThick',
                'differenceThick{}ContourValues'.format(hemisphere))
            differenceColormap = plt.get_cmap(config.get(
                'regriddedSeaIceConcThick', 'differenceColormap'))
            differenceColormapIndices = config.getExpression(
                'regriddedSeaIceConcThick', 'differenceColormapIndices')
            differenceColormap = cols.ListedColormap(
                differenceColormap(differenceColormapIndices),
                "differenceColormap")

            referenceLongitude = config.getfloat(
                'regriddedSeaIceConcThick',
                'referenceLongitude{}'.format(hemisphere))
            minimumLatitude = config.getfloat(
                'regriddedSeaIceConcThick',
                'minimumLatitude{}'.format(hemisphere))

            obsFileNames = buildConfigFullPath(
                config, 'seaIceObservations',
                'thickness{}_{}'.format(hemisphere, months))
            if not os.path.isfile(obsFileNames):
                raise SystemExit("Obs file {} not found. Exiting...".format(
                    obsFileNames))

            ncFile = netCDF4.Dataset(obsFileNames, mode='r')
            obsIceThickness = ncFile.variables["HI"][:].T
            ncFile.close()

            # Mask thickness fields
            iceThickness[iceThickness == 0] = ma.masked
            obsIceThickness = ma.masked_values(obsIceThickness, 0)
            if hemisphere == 'NH':
                # Obs thickness should be nan above 86 (ICESat data)
                obsIceThickness[latTarg > 86] = ma.masked
                plotProjection = 'npstere'
            else:
                plotProjection = 'spstere'

            difference = iceThickness - obsIceThickness

            title = "{} ({}, years {:04d}-{:04d})".format(subtitle, months,
                                                          startYear, endYear)
            fileout = "{}/icethick{}_{}_{}_years{:04d}-{:04d}.png".format(
                plotsDirectory, hemisphere, mainRunName, months, startYear,
                endYear)
            plot_polar_comparison(
                config,
                lonTarg,
                latTarg,
                iceThickness,
                obsIceThickness,
                difference,
                resultColormap,
                resultThickContourValues,
                differenceColormap,
                differenceThickContourValues,
                title=title,
                fileout=fileout,
                plotProjection=plotProjection,
                latmin=minimumLatitude,
                lon0=referenceLongitude,
                modelTitle=mainRunName,
                obsTitle="Observations (ICESat)",
                diffTitle="Model-Observations",
                cbarlabel="m")


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
