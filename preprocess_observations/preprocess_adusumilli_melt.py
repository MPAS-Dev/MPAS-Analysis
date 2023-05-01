#!/usr/bin/env python
import os
import shutil
import subprocess

import h5py
import numpy as np
import pyproj
import xarray as xr
from geometric_features import GeometricFeatures
from geometric_features.aggregation import get_aggregator_by_name
from mpas_tools.cime.constants import constants as cime_constants

from mpas_analysis.shared.io.download import download_files
from mpas_analysis.shared.constants import constants as mpas_constants


def download_adusumilli(out_filename):
    """
    Remap the Adusumilli et al. (2020) melt rates at 1 km resolution to an MPAS
    mesh

    Parameters
    ----------
    out_filename : str
        The original Adusumilli et al. (2020) melt rates
    """

    if os.path.exists(out_filename):
        return

    download_files(fileList=['_3_1.h5'],
                   urlBase='http://library.ucsd.edu/dc/object/bb0448974g/',
                   outDir='.')
    shutil.move('_3_1.h5', out_filename)


def adusumilli_hdf5_to_netcdf(in_filename, out_filename):
    """
    Remap the Adusumilli et al. (2020) melt rates at 1 km resolution to an MPAS
    mesh

    Parameters
    ----------
    in_filename : str
        The original Adusumilli et al. (2020) melt rates

    out_filename : str
        The Adusumilli et al. (2020) melt rates in NetCDF format
    """
    if os.path.exists(out_filename):
        return

    print(f'Reading {in_filename}...')
    h5_data = h5py.File(in_filename, 'r')

    x = np.array(h5_data['/x'])[:, 0]
    y = np.array(h5_data['/y'])[:, 0]
    melt_rate = np.array(h5_data['/w_b'])
    melt_rate_uncertainty = np.array(h5_data['/w_b_uncert'])
    print('done.')

    print('Creating an xarray dataset...')
    ds = xr.Dataset()

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')
    latlon_projection = pyproj.Proj(proj='latlong', datum='WGS84')

    print('Computing lat/lon...')
    x_2d, y_2d = np.meshgrid(x, y)
    transformer = pyproj.Transformer.from_proj(projection,
                                               latlon_projection)
    lon, lat = transformer.transform(x_2d, y_2d)
    print('done.')

    # Adusumilli et al. (2020) ice density (caption of Fig. 1 and Methods
    # section)
    rho_ice = 917.
    rho_fw = cime_constants['SHR_CONST_RHOFW']
    ice_to_fw_equiv = rho_ice / rho_fw

    ds['x'] = (('x',), x)
    ds['y'] = (('y',), y)
    ds['lon'] = (('y', 'x'), lon)
    ds.lon.attrs['units'] = 'degrees'
    ds['lat'] = (('y', 'x'), lat)
    ds.lat.attrs['units'] = 'degrees'
    ds['meltRate'] = (('y', 'x'), ice_to_fw_equiv * melt_rate)
    ds.meltRate.attrs['units'] = 'm/yr of freshwater'
    ds['meltRateUncertainty'] = (('y', 'x'),
                                 ice_to_fw_equiv * melt_rate_uncertainty)
    ds.meltRateUncertainty.attrs['units'] = 'm/yr of freshwater'
    print('Writing the dataset...')
    ds.to_netcdf(out_filename)
    print('done.')


def compute_regional_means(in_filename, out_prefix, chunk_size=1000,
                           core_count=128, multiprocessing_method='spawn'):
    """
    Remap the Adusumilli et al. (2020) melt rates at 1 km resolution to an MPAS
    mesh

    Parameters
    ----------
    in_filename : str
        The Adusumilli et al. (2020) melt rates in NetCDF format

    out_prefix : str
        A prefix for the file to contain the Adusumilli et al. (2020) mean melt
        rates and melt fluxes aggregated over ice-shelf regions

    chunk_size : int, optional
        The number of grid points that are processed in one operation when
        creating region masks

    core_count : int, optional
        The number of processes to use to compute masks.

    multiprocessing_method : str, optional
        The multiprocessing method use for python mask creation ('fork',
        'spawn' or 'forkserver')
    """
    mask_filename = 'ice_shelf_masks.nc'
    region_group = 'Ice Shelves'
    aggregation_function, prefix, date = \
        get_aggregator_by_name(region_group)
    out_filename = f'{out_prefix}_{prefix}{date}.nc'
    if not os.path.exists(mask_filename):
        geojson_filename = f'{prefix}{date}.geojson'
        gf = GeometricFeatures()
        fc = aggregation_function(gf)
        fc.to_geojson(geojson_filename)

        args = ['compute_projection_region_masks',
                '-i', in_filename,
                '--lon', 'lon',
                '--lat', 'lat',
                '-g', geojson_filename,
                '-o', mask_filename,
                '--chunk_size', f'{chunk_size}',
                '--process_count', f'{core_count}',
                '--multiprocessing_method', f'{multiprocessing_method}',
                '--show_progress']
        subprocess.run(args=args, check=True)

    ds = xr.open_dataset(in_filename)
    cell_area = ((ds.x.values[1] - ds.x.values[0]) *
                 (ds.y.values[1] - ds.y.values[0]))

    rho_fw = cime_constants['SHR_CONST_RHOFW']
    kg_per_gt = mpas_constants.kg_per_GT
    gt_per_m3 = rho_fw / kg_per_gt

    ds_out = xr.Dataset()
    ds_mask = xr.open_dataset(mask_filename)
    ds_out['regionNames'] = ds_mask.regionNames
    mean_melt_rate = []
    total_melt_flux = []
    area_melt = []
    melt_rate_uncert = []
    melt_flux_uncert = []
    area_uncert = []
    for region_index, region_name in enumerate(ds_mask.regionNames.values):
        region_mask = ds_mask.regionMasks.isel(nRegions=region_index) == 1

        melt_mean, melt_area = _compute_mean_and_area(
            ds.meltRate.values, region_mask, cell_area)
        melt_total = gt_per_m3 * melt_area * melt_mean
        mean_melt_rate.append(melt_mean)
        area_melt.append(melt_area)
        total_melt_flux.append(melt_total)
        uncert_sqr = ds.meltRateUncertainty.values**2
        uncert_mean_sqr, uncert_area = _compute_mean_and_area(
            uncert_sqr, region_mask, cell_area)
        uncert_rms = np.sqrt(uncert_mean_sqr)
        area_uncert.append(uncert_area)
        uncert_total = gt_per_m3 * uncert_area * uncert_rms
        melt_rate_uncert.append(uncert_rms)
        melt_flux_uncert.append(uncert_total)

    ds_out['meanMeltRate'] = (('nRegions',), mean_melt_rate)
    ds_out.meanMeltRate.attrs['units'] = 'm/yr of freshwater'
    ds_out['meltRateUncertainty'] = (('nRegions',), melt_rate_uncert)
    ds_out.meltRateUncertainty.attrs['units'] = 'm/yr of freshwater'
    ds_out['totalMeltFlux'] = (('nRegions',), total_melt_flux)
    ds_out.totalMeltFlux.attrs['units'] = 'GT/yr'
    ds_out['meltFluxUncertainty'] = (('nRegions',), melt_flux_uncert)
    ds_out.meltFluxUncertainty.attrs['units'] = 'GT/yr'
    ds_out['meltArea'] = (('nRegions',), area_melt)
    ds_out.meltArea.attrs['units'] = 'm^2'
    ds_out['uncertaintyArea'] = (('nRegions',), area_uncert)
    ds_out.uncertaintyArea.attrs['units'] = 'm^2'
    ds_out.to_netcdf(out_filename)


def _compute_mean_and_area(field, region_mask, cell_area):
    valid_melt = np.isfinite(field)
    mask = np.logical_and(region_mask, valid_melt)
    valid_count = np.count_nonzero(mask)
    area = cell_area * valid_count
    field_mean = np.sum(field[mask] * cell_area) / area
    return field_mean, area


def main():
    prefix = 'Adusumilli_2020_iceshelf_melt_rates_2010-2018_v0'
    hdf5_filename = f'{prefix}.h5'
    netcdf_filename = f'{prefix}.nc'

    download_adusumilli(hdf5_filename)
    adusumilli_hdf5_to_netcdf(hdf5_filename, netcdf_filename)
    compute_regional_means(netcdf_filename, prefix)


if __name__ == '__main__':
    main()
