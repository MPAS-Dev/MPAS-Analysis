"""
Computation and plotting of model/observation ocean transects.

Authors
-------
Jeremy Fyke, Xylar Asay-Davis

Last Modified
-------------
04/23/2017
"""

import numpy as np
import xarray as xr
import os
from pyproj import Geod, Proj
from scipy.interpolate import LinearNDInterpolator, griddata
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.ticker as ticker

from ..shared.constants import constants

from ..shared.io.utility import build_config_full_path

from ..shared.generalized_reader.generalized_reader \
    import open_multifile_dataset

from ..shared.timekeeping.utility import get_simulation_start_time

from ..shared.climatology.climatology import get_mpas_climatology_file_names, \
    cache_climatologies

from ..shared.analysis_task import setup_task

from ..shared.mpas_xarray import mpas_xarray


def antarctic_transects(config, fieldName):
    """
    Plots transects of ACME/MPAS output around the Antarctic continental shelf
    and compares them with reanalysis data

    Parameters
    ----------
    config :  instance of MpasAnalysisConfigParser
        Contains configuration options

    fieldName : {'temperature', 'salinity'}
        The name of a field to be analyized

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/22/2017
    """

    # perform common setup for the task
    namelist, runStreams, historyStreams, calendar, namelistMap, streamMap, \
        variableMap, plotsDirectory = setup_task(config, componentName='ocean')

    simulationStartTime = get_simulation_start_time(runStreams)

    # get a list of timeSeriesStats output files from the streams file,
    # reading only those that are between the start and end dates
    startDate = config.get('climatology', 'startDate')
    endDate = config.get('climatology', 'endDate')
    streamName = historyStreams.find_stream(streamMap['timeSeriesStats'])
    inputFiles = historyStreams.readpath(streamName, startDate=startDate,
                                         endDate=endDate, calendar=calendar)
    print '\n  Reading files:\n' \
          '    {} through\n    {}'.format(
              os.path.basename(inputFiles[0]),
              os.path.basename(inputFiles[-1]))

    observationsDirectory = build_config_full_path(
        config, 'oceanObservations', 'antarcticTransectsSubdirectory')

    try:
        restartFileName = runStreams.readpath('restart')[0]
    except ValueError:
        raise IOError('No MPAS-O restart file found: need at least one '
                      'restart file for antarctic_transects calculation')

    outputTimes = config.getExpression('antarcticTransects', 'comparisonTimes')

    # load MPAS data
    print("   Loading MPAS data...")

    dsMpasMesh = _load_mpas_mesh_data(restartFileName)

    varList = [fieldName, 'layerThickness', 'zMid']
    dsMpasData = open_multifile_dataset(
        fileNames=inputFiles, calendar=calendar, config=config,
        simulationStartTime=simulationStartTime,  timeVariableName='Time',
        variableList=varList, variableMap=variableMap, startDate=startDate,
        endDate=endDate)

    # mask out invalid values in the field (below the bathymetry)
    if 'zMid' in dsMpasData:
        # we don't need layer thickness if we already have zMid
        dsMpasData.drop('layerThickness')

    dsObs = _load_observations(config, fieldName, observationsDirectory)

    transectDepth, transectDistance, transectLonLats, nTransectPoints, \
        transectResoltuion = _define_transects(config)

    for monthsName in outputTimes:
        climatologyFileName, climatologyPrefix, regriddedFileName = \
            get_mpas_climatology_file_names(config, fieldName, monthsName)
        mpasClimatology = cache_climatologies(
            dsMpasData, constants.monthDictionary[monthsName], config,
            climatologyPrefix, calendar, printProgress=True)

        mpasClimatology[fieldName] = \
            mpasClimatology[fieldName].where(dsMpasMesh.cellMask)

        if 'zMid' not in mpasClimatology:
            mpasClimatology.load()
            mpasClimatology.close()
            # compute zMid from bottomDepth and layerThickness
            zMid = -mpasClimatology.layerThickness.cumsum(dim='nVertLevels') \
                + 0.5*mpasClimatology.layerThickness
            # match up the bottom with -bottomDepth
            zMid -= dsMpasMesh.bottomDepth + zMid.isel(nVertLevels=-1)
            mpasClimatology['zMid'] = zMid

            mpasClimatology.to_netcdf(climatologyFileName)

        _plot_transects(config, transectDepth, transectDistance,
                        transectLonLats, nTransectPoints, transectResoltuion,
                        dsMpasMesh, mpasClimatology, dsObs, fieldName,
                        plotsDirectory, monthsName)


def _load_mpas_mesh_data(restartFileName):
    """
    Load latitude, longitude and bottom depth, and compute cell mask from MPAS
    mesh file. Results are returned as xarray.DataArray objects

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/23/2017
    """

    dsMesh = xr.open_dataset(restartFileName)
    nCells = dsMesh.dims['nCells']
    nVertLevels = dsMesh.dims['nVertLevels']
    maxLevelCell = dsMesh.maxLevelCell.values - 1  # 1-based to 0-based index
    cellMask = np.zeros((nCells, nVertLevels), bool)
    for zIndex in range(nVertLevels):
        cellMask[zIndex <= maxLevelCell, zIndex] = True

    dsMesh = mpas_xarray.subset_variables(dsMesh, ['latCell', 'lonCell',
                                                   'bottomDepth'])

    dsMesh['cellMask'] = xr.DataArray(cellMask, dims=('nCells', 'nVertLevels'))
    return dsMesh


def _load_observations(config, fieldName, observationsDirectory):
    """
    Load and time average SOSE reanalysis data as an xarray data set

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/25/2017
    """

    years = '2005-2010'
    if fieldName == 'temperature':
        description = 'SOSE {} annual mean potential ' \
                      'temperature'.format(years)
        units = '$^\circ$C'
    elif fieldName == 'salinity':
        description = 'SOSE {} annual mean ' \
                      'salinity'.format(years)
        units = 'PSU'
    else:
        raise ValueError('Unexpected fieldName {}'.format(fieldName))

    inFileName = '{}/{}_SOSE_0.167x0.167_years{}.nc'.format(
        observationsDirectory, fieldName, years)

    ds = xr.open_dataset(inFileName)

    Lon, Lat = np.meshgrid(np.deg2rad(ds.lon.values),
                           np.deg2rad(ds.lat.values), indexing='ij')
    Lon = np.ravel(Lon)
    Lat = np.ravel(Lat)
    nCells = len(Lon)
    depth = ds.depth.values
    nDepth = len(depth)
    depth = np.tile(np.reshape(depth, (1, nDepth)), (nCells, 1))

    field = np.reshape(ds[fieldName].values, (nCells, nDepth))

    landMask = np.isnan(field[:, 0])

    field = field[-landMask, :]  # ... and remove these cells from SOSE
    depth = depth[-landMask, :]
    Lon = Lon[-landMask]
    Lat = Lat[-landMask]
    nCells = np.size(Lon)  # revise # of SOSE cells to just ocean

    dictonary = {'dims': ['nCells', 'nDepth'],
                 'coords': {'lon':
                            {'dims': ('nCells'),
                             'data': Lon,
                             'attrs': {'units': 'radians'}},
                            'lat':
                            {'dims': ('nCells'),
                             'data': Lat,
                             'attrs': {'units': 'radians'}},
                            'depth':
                            {'dims': ('nCells', 'nDepth'),
                             'data': depth,
                             'attrs': {'units': 'm'}}},
                 'data_vars': {fieldName:
                               {'dims': ('nCells', 'nDepth'),
                                'data': field,
                                'attrs': {'units': units,
                                          'description': description}}}}
    ds = xr.Dataset.from_dict(dictonary)

    # probably no need to cache this file, since the changes are small from
    # the observationos already stored and the data would be hard to interpret
    # ds.to_netcdf('{}_SOSE_DEBUG.nc'.format(fieldName))
    return ds


def _define_transects(config):
    """
    Read bounds from config file and construct a lists of transect properties

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/23/2017
    """

    # define transects here.  This could be done elsewhere.
    print("    Making transects...")

    sectionName = 'antarcticTransects'

    geod = Geod(ellps='sphere')
    # convert resoltion from km to m
    transectResolution = 1e3*config.getfloat(sectionName, 'transectResolution')
    transectDepth = config.getExpression(sectionName, 'transectDepth',
                                         usenumpyfunc=True)
    transectLat = config.getExpression(sectionName, 'transectLat',
                                       usenumpyfunc=True)
    transectLon = config.getExpression(sectionName, 'transectLon',
                                       usenumpyfunc=True)
    transectDistance = list()
    transectLonLats = list()
    nTransectPoints = list()
    for lon in transectLon:
        endPoints = np.deg2rad([transectLat[0], lon, transectLat[1], lon])

        # for code logic, start lon/lat must be <= end lon/lat
        if(endPoints[0] > endPoints[2]) or (endPoints[1] > endPoints[3]):
            raise ValueError("Error: start lon/lat > end lon/lat.")

        # Construct transect points in a great circle from start to end
        # lat/lon points
        distance = geod.inv(endPoints[1], endPoints[0],
                            endPoints[3], endPoints[2], radians=True)
        nPoints = np.floor(distance[2]/transectResolution).astype(int)
        transectPoints = geod.npts(endPoints[1], endPoints[0],
                                   endPoints[3], endPoints[2], nPoints,
                                   radians=True)

        # append transect info to lists of transect info
        nTransectPoints.append(nPoints)
        transectDistance.append(distance[2])
        transectLonLats.append(transectPoints)

    return transectDepth, transectDistance, transectLonLats, nTransectPoints, \
        transectResolution


def _generate_transect(transectLonLats, lon, lat, v, z, nCells, nDepth,
                       transectResolution, transectDepth, printDiags=False):
    """
    Generate a single MPAS or SOSE transect

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/23/2017
    """
    # put lon between -180 and 180
    lon = np.mod(lon + np.pi, 2*np.pi) - np.pi

    lonTransect, latTransect = zip(*transectLonLats)  # unpack lon,lat lists
    lonTransect = np.asarray(lonTransect)
    latTransect = np.asarray(latTransect)
    lats = latTransect[0]
    lons = lonTransect[0]
    late = latTransect[-1]
    lone = lonTransect[-1]
    nTransectPoints = np.size(lonTransect)
    # Vector with indices to transect points
    transectHasData = np.arange(nTransectPoints)

    if printDiags:
        print("Start/end transect lon: {}/{}".format(np.rad2deg(lons),
                                                     np.rad2deg(lone)))
        print("Start/end transect lat: {}/{}".format(np.rad2deg(lats),
                                                     np.rad2deg(late)))
        print("Min/max lon: {}/{}".format(np.rad2deg(np.min(lon)),
                                          np.rad2deg(np.max(lon))))
        print("Min/max lat: {}/{}".format(np.rad2deg(np.min(lat)),
                                          np.rad2deg(np.max(lat))))

    interpBuffer = np.deg2rad(1.)  # 1 degree, in radians
    indices = np.where((lon >= lons - interpBuffer) &
                       (lon <= lone + interpBuffer) &
                       (lat >= lats - interpBuffer) &
                       (lat <= late + interpBuffer))
    if printDiags:
        print("Total # of interpolant points: " + str(np.size(indices)))

    # Get minimu/maximum range of ocean values
    latMasked = lat[indices]
    lonMasked = lon[indices]

    latmin = np.min(latMasked)
    latmax = np.max(latMasked)
    latmid = np.mean([latmin, latmax])
    lonmin = np.min(lonMasked)
    lonmax = np.max(lonMasked)
    lonmid = np.mean([lonmin, lonmax])
    # Make regional projection here
    pgrid = Proj(proj='stere', lat_0=np.rad2deg(latmid),
                 lon_0=np.rad2deg(lonmid), ellps="sphere")

    if printDiags:
        print("Mid near-transect lon/lat: {}/{}".format(np.rad2deg(latmid),
                                                        np.rad2deg(lonmid)))
        print("Min/max near-transect lon: {}/{}".format(np.rad2deg(lonmin),
                                                        np.rad2deg(lonmax)))
        print("Min/max near-transect lat: {}/{}".format(np.rad2deg(latmin),
                                                        np.rad2deg(latmax)))
        print("Total # of interpolant points prior to coastal trimming: " +
              str(np.size(latTransect)))

    # Find transect points that fall within min/max lat/lons
    transectIndices = np.where((latTransect > latmin) &
                               (latTransect < latmax) &
                               (lonTransect > lonmin) &
                               (lonTransect < lonmax))
    lonTransect = lonTransect[transectIndices]
    latTransect = latTransect[transectIndices]
    # Indices with no-data points removed
    transectHasData = transectHasData[transectIndices]
    # Revise number of transect points
    nTransectPoints = np.size(lonTransect)

    if printDiags:
        print("Total # of interpolant points after coastal trimming: " +
              str(nTransectPoints))
    x, y = pgrid(np.rad2deg(lonMasked), np.rad2deg(latMasked))

    points = np.column_stack((x, y))
    xTransect, yTransect = pgrid(np.rad2deg(lonTransect),
                                 np.rad2deg(latTransect))
    transect = np.zeros((nDepth, nTransectPoints))
    depth = np.zeros((nDepth, nTransectPoints))
    for d in range(nDepth):  # iterate over depths
        VariableInterpolator = LinearNDInterpolator(points,
                                                    np.squeeze(v[indices, d]))
        DepthInterpolator = LinearNDInterpolator(points,
                                                 np.squeeze(z[indices, d]))
        transect[d, :] = VariableInterpolator(xTransect, yTransect)
        depth[d, :] = DepthInterpolator(xTransect, yTransect)
    mask = np.any(np.isnan(depth), axis=0)
    depth = depth[:, -mask]
    transect = transect[:, -mask]
    transectHasData = transectHasData[-mask]  # Indices with no data removed
    nTransectPoints = transect.shape[1]
    if printDiags:
        print("Total # of interpolant points after NaN trimming: " +
              str(nTransectPoints))

    xi = np.arange(nTransectPoints)*transectResolution  # along-transect points
    yi = transectDepth  # depth
    grid_x, grid_y = np.meshgrid(xi, yi)
    x = np.tile(xi, (nDepth, 1))
    transect_gridded = griddata((np.ravel(x), np.ravel(depth)),
                                np.ravel(transect), (grid_x, grid_y),
                                method='linear')

    return transect_gridded, xi, yi, transectHasData


def _plot_transects(config, transectDepth, transectDistance, transectLonLats,
                    nTransectPoints, transectResolution, dsMpasMesh,
                    dsMpasData, dsSose, fieldName, plotsDirectory, monthsName,
                    printDiags=False):
    """
    Plot images comparing MPAS to SOSE transects

    Authors
    -------
    Jeremy Fyke, Xylar Asay-Davis

    Last Modified
    -------------
    04/23/2017
    """
    if fieldName == 'temperature':
        cbarLabel = 'Temperature ($^\circ$C)'
        filePrefix = '{}/antarcticTransectTemp_{}'.format(plotsDirectory,
                                                          monthsName)

    elif fieldName == 'salinity':
        cbarLabel = 'Salinity (PSU)'
        filePrefix = '{}/antarcticTransectSalin_{}'.format(plotsDirectory,
                                                           monthsName)
    else:
        raise ValueError('Unexpected fieldName {}'.format(fieldName))

    resultLimits = config.getExpression(
        'antarcticTransects', '{}ResultLimits'.format(fieldName))
    differenceLimits = config.getExpression(
        'antarcticTransects', '{}DifferenceLimits'.format(fieldName))

    nTransects = len(transectDistance)

    # for each transect location, generate some transects (call to
    # _generate_transect)
    print("   Calling transect generator...")
    for t in range(nTransects):
        print("")
        lon = np.rad2deg(transectLonLats[t][0][0])
        print("    ***Generating transect {} at lon. {} ***".format(t, lon))
        if printDiags:
            print("     SOSE transect:")
        SOSEtransect, SOSEx, SOSEy, SOSEtransectHasData = \
            _generate_transect(transectLonLats[t], dsSose.lon.values,
                               dsSose.lat.values, dsSose[fieldName].values,
                               dsSose.depth.values, dsSose.dims['nCells'],
                               dsSose.dims['nDepth'], transectResolution,
                               transectDepth)

        if printDiags:
            print("     MPAS transect:")
        MPAStransect, MPASx, MPASy, MPAStransectHasData = \
            _generate_transect(transectLonLats[t], dsMpasMesh.lonCell.values,
                               dsMpasMesh.latCell.values,
                               dsMpasData[fieldName].values,
                               dsMpasData.zMid.values,
                               dsMpasData.dims['nCells'],
                               dsMpasData.dims['nVertLevels'],
                               transectResolution, transectDepth)

        # calculate difference between transects.  Since coastlines differ
        # between input datasets, *transectHasData is used to identify where
        # data exists for both datasets, and differences can be calculated
        diff = np.zeros((len(transectDepth), nTransectPoints[t]))
        diff[:, :] = np.nan
        minIndex = 0
        minIndexSet = False
        for n in range(nTransectPoints[t]):
            iSOSE = np.where(SOSEtransectHasData == n)
            iMPAS = np.where(MPAStransectHasData == n)
            if np.asarray(iSOSE).size + np.asarray(iMPAS).size == 2:
                # if both transects have same index
                if not minIndexSet:
                    minIndex = n
                    minIndexSet = True
                diff[:, n] = np.squeeze(MPAStransect[:, iMPAS]) - \
                    np.squeeze(SOSEtransect[:, iSOSE])
        diff = diff[:, minIndex:-1]
        dm = np.shape(diff)
        ntransectPoints_trimmed = dm[1]

        # do some plotting.
        levels = np.linspace(resultLimits[0], resultLimits[1], 50)
        fig = plt.figure()
        fig.set_size_inches(11, 8)
        ax1 = plt.subplot2grid((10, 3), (3, 0), rowspan=2, colspan=3)
        cs = ax1.contourf(SOSEx, SOSEy, SOSEtransect, levels,
                          cmap=plt.get_cmap("coolwarm"), extend="both")
        ticks_format = ticker.FuncFormatter(
            lambda x, pos: '{0:g}'.format(x/1000.))  # make axis km scale, P1
        ax1.xaxis.set_major_formatter(ticks_format)  # make axis km scale, P2
        ax1.yaxis.set_major_formatter(ticks_format)  # make axis km scale, P3
        ax1.set_xlabel('Distance (S to N, km)')
        ax1.set_ylabel('Depth (km)')
        ax1.set_title('SOSE transect')
        cbar = fig.colorbar(cs, ax=ax1, ticks=[resultLimits[0], 0,
                                               resultLimits[1]])
        cbar.ax.set_ylabel(cbarLabel, rotation=270)

        ax2 = plt.subplot2grid((10, 3), (0, 0), rowspan=2, colspan=3)
        cs = ax2.contourf(MPASx, MPASy, MPAStransect, levels,
                          cmap=plt.get_cmap("coolwarm"), extend="both")
        ax2.xaxis.set_major_formatter(ticks_format)  # make axis km scale, P2
        ax2.yaxis.set_major_formatter(ticks_format)  # make axis km scale, P3
        ax2.set_xlabel('Distance (S to N, km)')
        ax2.set_ylabel('Depth (km)')
        ax2.set_title('MPAS transect at {} longitude'.format(lon))
        cbar = fig.colorbar(cs, ax=ax2, ticks=[resultLimits[0], 0,
                                               resultLimits[1]])
        cbar.ax.set_ylabel(cbarLabel, rotation=270)

        levels = np.linspace(differenceLimits[0], differenceLimits[1], 50)
        xi = np.arange(ntransectPoints_trimmed) * transectResolution
        yi = transectDepth
        ax3 = plt.subplot2grid((10, 3), (6, 0), rowspan=2, colspan=3)
        cs = ax3.contourf(xi, yi, diff, levels, cmap=plt.get_cmap("coolwarm"),
                          extend="both")
        ax3.xaxis.set_major_formatter(ticks_format)  # make axis km scale, P2
        ax3.yaxis.set_major_formatter(ticks_format)  # make axis km scale, P3
        cbar = fig.colorbar(cs, ax=ax3, ticks=[differenceLimits[0], 0,
                                               differenceLimits[1]])
        cbar.ax.set_ylabel(cbarLabel, rotation=270)
        ax3.set_xlabel('Distance (S to N, km)')
        ax3.set_ylabel('Depth (km)')
        ax3.set_title('Difference (MPAS - SOSE)')

        # unpack lon,lat lists
        lonTransect, latTransect = zip(*transectLonLats[t])
        lonTransect = np.asarray(lonTransect) + np.pi  # convert to 0->2pi
        latTransect = np.asarray(latTransect)
        ax4 = plt.subplot2grid((10, 3), (9, 0), colspan=3)
        m = Basemap(llcrnrlon=0, llcrnrlat=-85, urcrnrlon=360, urcrnrlat=-50,
                    resolution='l', ax=ax4)
        m.drawcoastlines(linewidth=0.25)
        m.fillcontinents(color='black', lake_color='aqua')
        # draw the edge of the map projection region (the projection limb)
        m.drawmapboundary(fill_color='aqua')
        # draw lat/lon grid lines every 30 degrees.
        m.drawmeridians(np.arange(0, 360, 30))
        m.drawparallels(np.arange(-90, 90, 30))
        x, y = m([np.rad2deg(lonTransect[0]), np.rad2deg(lonTransect[-1])],
                 [np.rad2deg(latTransect[0]), np.rad2deg(latTransect[-1])])
        m.plot(x, y, 'r')

        plt.savefig("{}_transect_{:03}_lon_{}.png".format(filePrefix, t, lon),
                    bbox_inches='tight')
        plt.close(fig)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
