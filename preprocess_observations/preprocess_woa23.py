#!/usr/bin/env python
import os
from datetime import datetime

import gsw
import numpy as np
import xarray as xr
from mpas_analysis.shared.climatology import compute_climatology
from mpas_analysis.shared.constants import constants
from mpas_analysis.shared.io.download import download_files
from mpas_analysis.shared.io import write_netcdf
from mpas_tools.cime.constants import constants as cime_constants


def download_woa():
    """
    Download 0.25 degree WOA23 annual and monthly files
    """
    base_url = \
        'https://www.ncei.noaa.gov/thredds-ocean/fileServer/woa23/DATA'

    for field in ['temperature', 'salinity']:
        for index in range(13):
            woa_filename = f'woa23_decav91C0_{field[0]}{index:02d}_04.nc'
            woa_url = f'{base_url}/{field}/netcdf/decav91C0/0.25/'

            if os.path.exists(woa_filename):
                continue

            download_files(fileList=[woa_filename],
                           urlBase=woa_url,
                           outDir='.')


def combine():
    """
    Run this step of the test case
    """
    now = datetime.now()
    datestring = now.strftime("%Y%m%d")
    out_filename = f'woa23_decav_04_pt_s_mon_ann.{datestring}.nc'
    if os.path.exists(out_filename):
        return

    ds_ann = xr.open_dataset('woa23_decav91C0_t00_04.nc', decode_times=False)

    ds_out = xr.Dataset()

    for var in ['lon', 'lat', 'depth']:
        ds_out[var] = ds_ann[var]
        ds_out[f'{var}_bnds'] = ds_ann[f'{var}_bnds']

    for field in ['temperature', 'salinity']:
        var = f'{field[0]}_an'
        ann_filename = f'woa23_decav91C0_{field[0]}00_04.nc'
        ds_ann = xr.open_dataset(
            ann_filename, decode_times=False).isel(time=0).drop_vars('time')

        time_slices = []
        for month in range(1, 13):
            depth_slices = []
            ds_month = xr.open_dataset(
                f'woa23_decav91C0_{field[0]}{month:02d}_04.nc',
                decode_times=False).isel(time=0).drop_vars('time')
            for depth_index in range(ds_ann.sizes['depth']):
                if depth_index < ds_month.sizes['depth']:
                    ds = ds_month
                else:
                    ds = ds_ann
                depth_slices.append(ds[var].isel(depth=depth_index))

            da_month = xr.concat(depth_slices, dim='depth')
            print(month, da_month)
            time_slices.append(da_month)

        ds_out[var] = xr.concat(time_slices, dim='time')
        ds_out[var] = ds_out[var].transpose('time', 'lat', 'lon', 'depth')
        ds_out[var].attrs = ds_ann[var].attrs

    ds_out['month'] = ('time', np.arange(1, 13))
    ds_out = _temp_to_pot_temp(ds_out)
    ds_out.to_netcdf(out_filename)


def compute_obs_ts_climatology(season):
    """
    season : str
        The season over which to compute the climatology
    """
    now = datetime.now()
    datestring = now.strftime("%Y%m%d")

    in_filename = f'woa23_decav_04_pt_s_mon_ann.{datestring}.nc'
    out_filename = f'woa23_{season}_decav_04_pt_s_z_vol.{datestring}.nc'

    if os.path.exists(out_filename):
        return

    print("\n computing T S climatogy for WOA23...")

    print(f'  Reading from {in_filename}...')
    ds = xr.open_dataset(in_filename)

    # compute volume from lat, lon, depth bounds
    print('  Computing volume...')
    lat_bnds = ds['lat_bnds']
    lon_bnds = ds['lon_bnds']
    depth_bnds = ds['depth_bnds']
    dlat = np.deg2rad(lat_bnds[:, 1] - lat_bnds[:, 0])
    dlon = np.deg2rad(lon_bnds[:, 1] - lon_bnds[:, 0])
    lat = np.deg2rad(ds['lat'])
    dz = depth_bnds[:, 1] - depth_bnds[:, 0]
    radius = cime_constants['SHR_CONST_REARTH']
    area = radius**2*np.cos(lat)*dlat*dlon
    ds['volume'] = dz*area

    ds = ds.rename({'time': 'Time'})
    ds = ds.set_coords('month')
    ds.coords['year'] = np.ones(ds.sizes['Time'], int)

    month_values = constants.monthDictionary[season]
    print('  Computing climatology...')
    ds = compute_climatology(ds, month_values, maskVaries=True)

    print('  Broadcasting z coordinate...')
    attrs = ds['depth'].attrs
    ds['depth'] = -ds['depth']
    ds['depth'].attrs = attrs
    ds['depth'].attrs['positive'] = 'up'

    _, _, z = xr.broadcast(ds['pt_an'], ds['s_an'], ds['depth'])

    ds['zBroadcast'] = z

    print(f'  writing {out_filename}...')
    write_netcdf.write_netcdf_with_fill(ds, out_filename)
    print('  Done!')


def main():
    download_woa()
    combine()
    for season in ['ANN', 'JFM', 'JAS']:
        compute_obs_ts_climatology(season)


def _temp_to_pot_temp(ds):
    dims = ds.t_an.dims

    slices = list()
    for depth_index in range(ds.sizes['depth']):
        temp_slice = ds.t_an.isel(depth=depth_index)
        in_situ_temp = temp_slice.values
        salin = ds.s_an.isel(depth=depth_index).values
        lat = ds.lat.broadcast_like(temp_slice).values
        lon = ds.lon.broadcast_like(temp_slice).values
        z = -ds.depth.isel(depth=depth_index).values
        pressure = gsw.p_from_z(z, lat)
        mask = np.isfinite(in_situ_temp)
        SA = gsw.SA_from_SP(salin[mask], pressure[mask], lon[mask],
                            lat[mask])
        pot_temp = np.nan * np.ones(in_situ_temp.shape)
        pot_temp[mask] = gsw.pt_from_t(SA, in_situ_temp[mask],
                                       pressure[mask], p_ref=0.)
        pot_temp_slice = xr.DataArray(data=pot_temp, dims=temp_slice.dims,
                                      attrs=temp_slice.attrs)

        slices.append(pot_temp_slice)

    ds['pt_an'] = xr.concat(slices, dim='depth').transpose(*dims)

    ds.pt_an.attrs['standard_name'] = \
        'sea_water_potential_temperature'
    ds.pt_an.attrs['long_name'] = \
        'Objectively analyzed mean fields for ' \
        'sea_water_potential_temperature at standard depth levels.'

    ds = ds.drop_vars('t_an')
    return ds


if __name__ == '__main__':
    main()
