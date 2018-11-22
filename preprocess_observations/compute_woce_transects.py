#!/usr/bin/env python

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import gsw
import xarray
import os
import numpy
import argparse
import zipfile
import glob
from datetime import datetime

from mpas_analysis.shared.io.utility import make_directories
from mpas_analysis.shared.io.download import download_files
from mpas_analysis.shared.io import write_netcdf


def process_transect(url, stations, outFileName, transectName, inDir, outDir):

    urlBase, file = url.rsplit('/', 1)
    download_files([file], urlBase=urlBase, outDir=inDir)

    with zipfile.ZipFile('{}/{}'.format(inDir, file), 'r') as f:
        f.extractall('{}/stations/'.format(inDir))

    inFileList = sorted(glob.glob('{}/stations/{}_*.nc'.format(
            inDir, transectName)))

    pressureValues = set()
    fileNames = {}
    validStations = set()

    for fileName in inFileList:
        path, file = os.path.split(fileName)
        _, station, cast, _ = file.split('_')
        station = int(station)
        cast = int(cast)
        if station not in stations:
            continue
        validStations.add(station)
        with xarray.open_dataset(fileName) as ds:
            for pressure in ds.pressure.values:
                pressureValues.add(pressure)

        if station not in fileNames:
            fileNames[station] = {}
        fileNames[station][cast] = fileName

    stations = validStations

    nStations = len(stations)
    nDepths = len(pressureValues)
    latitude = numpy.zeros(nStations)
    longitude = numpy.zeros(nStations)
    pressure = numpy.zeros((nStations, nDepths))
    temperature = numpy.zeros((nStations, nDepths))
    salinity = numpy.zeros((nStations, nDepths))
    nValues = numpy.zeros((nStations, nDepths), int)

    for stationIndex, station in enumerate(stations):

        castCount = len(fileNames[station])
        for cast in fileNames[station]:
            fileName = fileNames[station][cast]
            with xarray.open_dataset(fileName) as ds:
                nPoints = len(ds.pressure)
                nValues[stationIndex, 0:nPoints] += 1
                pressure[stationIndex, 0:nPoints] += ds.pressure
                temperature[stationIndex, 0:nPoints] += ds.temperature
                salinity[stationIndex, 0:nPoints] += ds.salinity

                longitude[stationIndex] += ds.longitude.values[0]
                latitude[stationIndex] += ds.latitude.values[0]

        longitude[stationIndex] /= castCount
        latitude[stationIndex] /= castCount

    # average over casts
    validMask = nValues > 0
    pressure = numpy.ma.masked_array(
            numpy.divide(pressure, nValues, where=validMask),
            mask=(nValues == 0))
    temperature = numpy.ma.masked_array(
            numpy.divide(temperature, nValues, where=validMask),
            mask=(nValues == 0))
    salinity = numpy.ma.masked_array(
            numpy.divide(salinity, nValues, where=validMask),
            mask=(nValues == 0))

    Lat = numpy.tile(latitude.reshape(nStations, 1), (1, nDepths))
    Lon = numpy.tile(longitude.reshape(nStations, 1), (1, nDepths))

    z = numpy.ma.masked_all(pressure.shape)
    z[validMask] = gsw.z_from_p(pressure[validMask], Lat[validMask])
    SA = gsw.SA_from_SP(salinity[validMask], pressure[validMask],
                        Lon[validMask], Lat[validMask])

    CT = gsw.CT_from_t(SA, temperature[validMask], pressure[validMask])

    potDensity = numpy.ma.masked_all(pressure.shape)
    potDensity[validMask] = gsw.rho(SA, CT, 0.)

    potTemp = numpy.ma.masked_all(pressure.shape)

    potTemp[validMask] = gsw.pt0_from_t(SA, temperature[validMask],
                                        pressure[validMask])

    longitude = xarray.DataArray.from_dict({'dims': ('nPoints',),
                                            'data': longitude,
                                            'attrs': {'long_name': 'longitude',
                                                      'units': 'degrees'}})
    latitude = xarray.DataArray.from_dict({'dims': ('nPoints',),
                                           'data': latitude,
                                           'attrs': {'long_name': 'latitude',
                                                     'units': 'degrees'}})

    pressure = xarray.DataArray.from_dict({'dims': ('nPoints',
                                                    'nz'),
                                           'data': pressure,
                                           'attrs':
                                               {'long_name': 'pressure',
                                                'units': 'dbar'}})

    z = xarray.DataArray.from_dict({'dims': ('nPoints', 'nz'),
                                    'data': z,
                                    'attrs':
                                        {'long_name': 'height',
                                         'units': 'm'}})

    potTemp = xarray.DataArray.from_dict({'dims': ('nPoints',
                                                   'nz'),
                                          'data': potTemp,
                                          'attrs':
                                              {'long_name':
                                                  'potential temperature',
                                               'units': 'deg C'}})

    salinity = xarray.DataArray.from_dict({'dims': ('nPoints',
                                                    'nz'),
                                           'data': salinity,
                                           'attrs':
                                               {'long_name': 'salinity',
                                                'units': 'PSU'}})

    potDensity = xarray.DataArray.from_dict(
            {'dims': ('nPoints', 'nz'),
             'data': potDensity,
             'attrs': {'long_name': 'potential denisty',
                       'units': 'kg m^{-3}'}})

    dsTransect = xarray.Dataset({'lon': longitude,
                                 'lat': latitude,
                                 'pressure': pressure,
                                 'z': z,
                                 'potentialTemperature': potTemp,
                                 'salinity': salinity,
                                 'potentialDensity': potDensity})

    write_netcdf(dsTransect, '{}/{}'.format(outDir, outFileName))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inDir", dest="inDir", required=True,
                        help="Directory where intermediate files used in "
                             "processing should be downloaded")
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where final preprocessed observation "
                             "are stored")
    args = parser.parse_args()

    make_directories(args.inDir)
    make_directories(args.outDir)

    date = datetime.now().strftime('%Y%m%d')

    process_transect(url='https://cchdo.ucsd.edu/data/4801/a21_nc_ctd.zip',
                     stations=numpy.arange(102, 118),
                     outFileName='WOCE_A21_Drake_Passage_{}.nc'.format(date),
                     transectName='a21',
                     inDir=args.inDir,
                     outDir=args.outDir)

    process_transect(url='https://cchdo.ucsd.edu/data/3654/a23_nc_ctd.zip',
                     stations=numpy.arange(3, 128),
                     outFileName='WOCE_A23_South_Atlantic_{}.nc'.format(date),
                     transectName='a23',
                     inDir=args.inDir,
                     outDir=args.outDir)

    process_transect(url='https://cchdo.ucsd.edu/data/4941/a12_nc_ctd.zip',
                     stations=numpy.arange(536, 607),
                     outFileName='WOCE_A12_Prime_Meridian_{}.nc'.format(date),
                     transectName='a12',
                     inDir=args.inDir,
                     outDir=args.outDir)
