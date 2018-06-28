#!/usr/bin/env python 
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

"""
A script for downloading and preprocessing BGC data sets for use in
MPAS-Analysis
"""
# Authors
# -------
# Riley X. Brady

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np
import xarray as xr
import os
import argparse
import glob
from mpas_analysis.shared.io.download import download_files


def process_landschuetzer(inDir, outDir):
    """
    Edit Landschuetzerv2016 raw netCDF output to conform with MPAS-Analysis
    infrastructure
    """
    # Authors
    # -------
    # Riley X. Brady
    print("Processing and saving CO2 flux data...")
    ds = xr.open_dataset(inDir + '/spco2_1982-2015_MPI_SOM-FFN_v2016.nc',
                         drop_variables='date')
    ds = ds['fgco2_smoothed']
    ds = ds.rename({'time': 'Time'})
    ds.name = 'CO2_gas_flux'
    ds.coords['month'] = ds['Time.month']
    ds.coords['year'] = ds['Time.year']
    outFile = ('CO2_gas_flux_1.0x1.0degree.nc')
    ds.to_dataset().to_netcdf(outDir + '/' + outFile)


def process_woa(variable, ref_char, inDir, outDir):
    """
    Edit WOA raw netCDF output to conform with MPAS-Analysis infrastructure
    """
    # Authors
    # -------
    # Riley X. Brady
    print("Processing and saving " + variable + "...")
    ds = xr.open_mfdataset(inDir + '/woa13_all_' + ref_char + '*.nc',
                           decode_times=False, concat_dim='time')
    ds = ds[ref_char + '_an'].isel(depth=0)
    ds.coords['month'] = ('time', np.arange(1, 13))
    ds.coords['year'] = ('time', np.ones(12))
    ds = ds.rename({'time': 'Time'})
    ds.name = variable
    del ds['Time']
    outFile = (variable + '_1.0x1.0degree.nc')
    ds.to_dataset().to_netcdf(outDir + '/' + outFile)


if __name__ == '__main__':
    """
    Process raw Landschuetzerv2016 and WOA observations for comparison to 
    MPAS output.
    """
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--inDir", dest="inDir", required=True,
                        help="Directory where intermediate files used in "
                             "processing should be downloaded")
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where final preprocessed observation "
                             "are stored")
    args = parser.parse_args()

    try:
        os.makedirs(args.inDir)
    except OSError:
        pass

    # + + + Landschuetzer Carbon Flux + + +
    urlBase = ('https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/' +
               'data/0-data/')
    download_files(['spco2_1982-2015_MPI_SOM-FFN_v2016.nc'],
                   urlBase=urlBase, outDir=args.inDir)
    # Skip processing and saving netCDF if it's already been done.
    if not os.path.isfile(args.outDir + '/' + 'CO2_gas_flux_1.0x1.0degree.nc'):
        process_landschuetzer(args.inDir, args.outDir)

    # + + +  World Ocean Atlas data + + +
    variables = ['NO3', 'PO4', 'SiO3', 'O2']
    # Character convention for WOA data
    ref_char = {'NO3': 'n', 'PO4': 'p', 'SiO3': 'i', 'O2': 'o'}
    # Naming convention for WOA data
    ref_name = {'NO3': 'nitrate', 'PO4': 'phosphate', 'SiO3': 'silicate',
                'O2': 'oxygen'}
    for v in variables:
        file_list = []
        for i in ["%.2d" % i for i in range(1, 13)]:
            file_list.append('woa13_all_' + ref_char[v] + i + '_01.nc')
        urlBase = ('https://data.nodc.noaa.gov/woa/WOA13/DATA/' + ref_name[v] +
                   '/netcdf/all/1.00/')
        download_files(file_list, urlBase=urlBase, outDir=args.inDir)
        # Skip processing and saving netCDF if it's already been done.
        if not os.path.isfile(args.outDir + '/' + v +
                              '_1.0x1.0degree.nc'):
            process_woa(v, ref_char[v], args.inDir, args.outDir)
     
