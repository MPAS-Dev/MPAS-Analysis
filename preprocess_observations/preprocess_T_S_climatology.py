#!/usr/bin/env python

import argparse
import xarray
import numpy

import dask
import multiprocessing
from multiprocessing.pool import ThreadPool

from mpas_analysis.shared.io import write_netcdf

from mpas_analysis.shared.constants import constants

from mpas_analysis.shared.climatology import compute_climatology


def compute_obs_ts_climatology(obs_name, obs_dict, base_dir, season):
    """
    Plots time-series output of properties in an ocean region.
    obs_name : str
        The name of the observational data set

    obs_dict : dict
        Information on the observational data sets

    base_dir : str
        The base directory where observations are found and the climatology
        should be written

    season : str
        The season over which to compute the climatology
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    componentName = 'ocean'

    fileName = obs_dict['outFileTemplate'].format(season)
    daskThreads = multiprocessing.cpu_count()

    base_dir = '{}/observations/Ocean'.format(base_dir)

    with dask.config.set(schedular='threads',
                         pool=ThreadPool(daskThreads)):

        print("\n computing T S climatogy for {}...".format(
            obs_name))

        chunk = {obs_dict['tVar']: 6}

        TVarName = obs_dict['TVar']
        SVarName = obs_dict['SVar']
        zVarName = obs_dict['zVar']
        lonVarName = obs_dict['lonVar']
        latVarName = obs_dict['latVar']
        volVarName = obs_dict['volVar']

        obsFileName = '{}/{}'.format(base_dir, obs_dict['TFileName'])
        print('  Reading from {}...'.format(obsFileName))
        ds = xarray.open_dataset(obsFileName, chunks=chunk)
        if obs_dict['SFileName'] != obs_dict['TFileName']:
            obsFileName = '{}/{}'.format(base_dir, obs_dict['SFileName'])
            print('  Reading from {}...'.format(obsFileName))
            dsS = xarray.open_dataset(obsFileName, chunks=chunk)
            ds[SVarName] = dsS[SVarName]

        if obs_dict['volFileName'] is None:
            # compute volume from lat, lon, depth bounds
            print('  Computing volume...'.format(obsFileName))
            latBndsName = ds[latVarName].attrs['bounds']
            lonBndsName = ds[lonVarName].attrs['bounds']
            zBndsName = ds[zVarName].attrs['bounds']
            latBnds = ds[latBndsName]
            lonBnds = ds[lonBndsName]
            zBnds = ds[zBndsName]
            dLat = numpy.deg2rad(latBnds[:, 1] - latBnds[:, 0])
            dLon = numpy.deg2rad(lonBnds[:, 1] - lonBnds[:, 0])
            lat = numpy.deg2rad(ds[latVarName])
            dz = zBnds[:, 1] - zBnds[:, 0]
            radius = 6378137.0
            area = radius**2*numpy.cos(lat)*dLat*dLon
            ds[volVarName] = dz*area

        elif obs_dict['volFileName'] != obs_dict['TFileName']:
            obsFileName = '{}/{}'.format(base_dir, obs_dict['volFileName'])
            print('  Reading from {}...'.format(obsFileName))
            dsVol = xarray.open_dataset(obsFileName)
            ds[volVarName] = dsVol[volVarName]

        chunk = {obs_dict['latVar']: 400,
                 obs_dict['lonVar']: 400}

        ds = ds.chunk(chunks=chunk)

        if obs_dict['tVar'] in ds.dims:
            if obs_dict['tVar'] != 'Time':
                if obs_dict['tVar'] == 'month':
                    ds = ds.rename({obs_dict['tVar']: 'Time'})
                    ds.coords['month'] = ds['Time']
                else:
                    ds = ds.rename({obs_dict['tVar']: 'Time'})
            if 'year' not in ds:
                ds.coords['year'] = numpy.ones(ds.sizes['Time'], int)

            monthValues = constants.monthDictionary[season]
            print('  Computing climatology...')
            ds = compute_climatology(ds, monthValues, maskVaries=True)

        print('  Broadcasting z coordinate...')
        if 'positive' in ds[zVarName].attrs and \
                ds[zVarName].attrs['positive'] == 'down':
            attrs = ds[zVarName].attrs
            ds[zVarName] = -ds[zVarName]
            ds[zVarName].attrs = attrs
            ds[zVarName].attrs['positive'] = 'up'

        T, S, z = xarray.broadcast(ds[TVarName], ds[SVarName],
                                   ds[zVarName])

        ds['zBroadcast'] = z

        print('  writing {}...'.format(fileName))
        fileName = '{}/{}'.format(base_dir, fileName)
        write_netcdf(ds, fileName)
        print('  Done!')


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-d", "--dir", dest="dir", required=True,
                        help="The base directory where diagnostics are stored")
    args = parser.parse_args()

    obs_dicts = {
        'SOSE': {
            'suffix': 'SOSE',
            'gridName': 'SouthernOcean_0.167x0.167degree',
            'gridFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                            'SouthernOcean_0.167x0.167degree_20180710.nc',
            'TFileName': 'SOSE/SOSE_2005-2010_monthly_pot_temp_'
                         'SouthernOcean_0.167x0.167degree_20180710.nc',
            'SFileName': 'SOSE/SOSE_2005-2010_monthly_salinity_'
                         'SouthernOcean_0.167x0.167degree_20180710.nc',
            'volFileName': 'SOSE/SOSE_volume_'
                           'SouthernOcean_0.167x0.167degree_20190815.nc',
            'outFileTemplate': 'SOSE/SOSE_{}_T_S_z_vol_'
                           'SouthernOcean_0.167x0.167degree_20200514.nc',
            'lonVar': 'lon',
            'latVar': 'lat',
            'TVar': 'theta',
            'SVar': 'salinity',
            'volVar': 'volume',
            'zVar': 'z',
            'tVar': 'Time'},
        'WOA18': {
            'suffix': 'WOA18',
            'gridName': 'Global_0.25x0.25degree',
            'gridFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'TFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'SFileName': 'WOA18/woa18_decav_04_TS_mon_20190829.nc',
            'volFileName': None,
            'outFileTemplate': 'WOA18/woa18_{}_T_S_z_vol_20200514.nc',
            'lonVar': 'lon',
            'latVar': 'lat',
            'TVar': 't_an',
            'SVar': 's_an',
            'volVar': 'volume',
            'zVar': 'depth',
            'tVar': 'month'}}

    seasons = ['ANN', 'JFM', 'JAS']

    for obs_name in obs_dicts:
        obs_dict = obs_dicts[obs_name]
        for season in seasons:
            compute_obs_ts_climatology(obs_name, obs_dict, args.dir, season)


if __name__ == '__main__':
    main()
