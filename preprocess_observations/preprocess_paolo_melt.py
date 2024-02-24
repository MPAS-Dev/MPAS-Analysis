#!/usr/bin/env python
import os
import shutil

import numpy as np
import pyproj
import xarray as xr
from mpas_tools.cime.constants import constants as cime_constants
from pyremap import Remapper, ProjectionGridDescriptor
from pyremap.polar import get_antarctic_stereographic_projection

from mpas_analysis.shared.io.download import download_files


def download_paolo(out_filename):
    """
    Remap the Paolo et al. (2023) melt rates at 1 km resolution to an MPAS
    mesh

    Parameters
    ----------
    out_filename : str
        The original Paolo et al. (2023) melt rates
    """

    if os.path.exists(out_filename):
        return

    download_files(fileList=['ANT_G1920V01_IceShelfMelt.nc'],
                   urlBase='https://its-live-data.s3.amazonaws.com/height_change/Antarctica/Floating',
                   outDir='.')
    shutil.move('ANT_G1920V01_IceShelfMelt.nc', out_filename)


def process_paolo(in_filename, out_filename):
    """
    Convert Paolo et al. (2023) melt rates to freshwater equivalent

    Parameters
    ----------
    in_filename : str
        The original Paolo et al. (2023) melt rates

    out_filename : str
        The Paolo et al. (2023) melt rates in freshwater equivalent
    """
    if os.path.exists(out_filename):
        return

    print(f'Reading {in_filename}...')

    with xr.open_dataset(in_filename) as ds_in:

        x = ds_in.x
        y = ds_in.y
        melt_rate = ds_in.melt_mean
        melt_rate_uncertainty = np.sqrt((ds_in.melt_err**2).mean(dim='time'))

    print('done.')

    print('Creating an xarray dataset...')
    ds = xr.Dataset()

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    latlon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

    print('Computing lat/lon...')
    x_2d, y_2d = np.meshgrid(x.values, y.values)
    transformer = pyproj.Transformer.from_proj(projection,
                                               latlon_projection)
    lon, lat = transformer.transform(x_2d, y_2d)
    print('done.')

    # Paolo et al. (2023) ice density
    rho_ice = 917.
    rho_fw = cime_constants['SHR_CONST_RHOFW']
    ice_to_fw_equiv = rho_ice / rho_fw

    ds['x'] = x
    ds['y'] = y
    ds['lon'] = (('y', 'x'), lon)
    ds.lon.attrs['units'] = 'degrees'
    ds['lat'] = (('y', 'x'), lat)
    ds.lat.attrs['units'] = 'degrees'
    ds['meltRate'] = -ice_to_fw_equiv * melt_rate
    ds.meltRate.attrs['units'] = 'm/yr of freshwater'
    ds['meltRateUncertainty'] = ice_to_fw_equiv * melt_rate_uncertainty
    ds.meltRateUncertainty.attrs['units'] = 'm/yr of freshwater'
    print('Writing the dataset...')
    ds.to_netcdf(out_filename)
    print('done.')


def remap_paolo(in_filename, out_prefix, date, task_count=128):
    """
    Remap Paolo et al. (2023) melt rates to comparison grids

    Parameters
    ----------
    in_filename : str
        The Paolo et al. (2023) melt rates in NetCDF format

    out_prefix : str
        A prefix for the file to contain the Paolo et al. (2023) melt
        rates and melt fluxes remapped to the comparison grid

    date : str
        A date string to append to the file name.

    task_count : int
        The number of MPI tasks to use to create the mapping file
    """
    ds = xr.open_dataset(in_filename)

    melt_attrs = ds.meltRate.attrs
    uncert_attrs = ds.meltRateUncertainty.attrs

    mask = ds.meltRate.notnull()
    ds['meltRate'] = ds.meltRate.where(mask, 0.)
    ds['meltMask'] = mask.astype(float)
    mask = ds.meltRateUncertainty.notnull()
    ds['meltRateUncertSqr'] = (ds.meltRateUncertainty**2).where(mask, 0.)
    ds['uncertMask'] = mask.astype(float)
    ds = ds.drop_vars(['lat', 'lon', 'meltRateUncertainty'])

    in_x = ds.x.values
    in_y = ds.y.values
    lx = np.abs(1e-3 * (in_x[-1] - in_x[0]))
    ly = np.abs(1e-3 * (in_y[-1] - in_y[0]))
    dx = np.abs(1e-3 * (in_x[1] - in_x[0]))

    in_grid_name = f'{lx:g}x{ly:g}km_{dx:g}km_Antarctic_stereo'

    in_projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 '
                                '+lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 '
                                '+ellps=WGS84')

    in_descriptor = ProjectionGridDescriptor.create(
        in_projection, in_x, in_y, in_grid_name)

    width = 6000.
    reses = [1., 4., 10.]

    for res in reses:
        x_max = 0.5 * width * 1e3
        nx = int(width / res) + 1
        out_x = np.linspace(-x_max, x_max, nx)

        out_grid_name = f'{width:g}x{width:g}km_{res:g}km_Antarctic_stereo'

        out_projection = get_antarctic_stereographic_projection()

        out_descriptor = ProjectionGridDescriptor.create(
            out_projection, out_x, out_x, out_grid_name)

        method = 'conserve'

        map_filename = f'map_{in_grid_name}_to_{out_grid_name}_{method}.nc'

        remapper = Remapper(in_descriptor, out_descriptor, map_filename)

        if not os.path.exists(map_filename):
            remapper.build_mapping_file(method=method, mpiTasks=task_count,
                                        esmf_parallel_exec='srun')

        ds_out = remapper.remap(ds)
        mask = ds_out.meltMask > 0.
        ds_out['meltRate'] = ds_out.meltRate.where(mask)
        ds_out.meltRate.attrs = melt_attrs
        mask = ds_out.uncertMask > 0.
        ds_out['meltRateUncertainty'] = \
            (np.sqrt(ds_out.meltRateUncertSqr)).where(mask)
        ds_out.meltRateUncertainty.attrs = uncert_attrs
        ds_out = ds_out.drop_vars(['meltRateUncertSqr'])
        ds_out.to_netcdf(f'{out_prefix}_{out_grid_name}.{date}.nc')


def main():
    prefix = 'Paolo_2023_iceshelf_melt_rates_1992-2017_v1.0'
    date = '20240220'

    orig_filename = 'Paolo_2023_ANT_G1920V01_IceShelfMelt.nc'
    processed_filename = f'{prefix}.{date}.nc'

    download_paolo(orig_filename)
    process_paolo(orig_filename, processed_filename)
    remap_paolo(processed_filename, prefix, date)


if __name__ == '__main__':
    main()
