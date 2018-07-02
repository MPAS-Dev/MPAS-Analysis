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
import tarfile

from mpas_analysis.shared.io.download import download_files
from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import LatLonGridDescriptor
from mpas_analysis.shared.climatology.comparison_descriptors \
    import get_comparison_descriptor
from mpas_analysis.configuration \
    import MpasAnalysisConfigParser

def process_SeaWIFS(inDir, outDir):
    """
    Load in all monthly climatology and process for MPAS-Analysis. 
    """
    # Authors
    # -------
    # Riley X. Brady
    if not os.path.isfile(outDir + '/Chl_SeaWIFS.nc'):
        print("Processing and saving SeaWIFS data...")
        ds = xr.open_mfdataset(inDir + '/' + 'S199*.nc', concat_dim='Time')
        ds = ds['chlor_a']
        ds['month'] = ('Time', range(1, 13))
        ds['year'] = ('Time', np.ones(12))
        ds.name = 'Chl'
        outFile = 'Chl_SeaWIFS.nc'
        ds.to_dataset().to_netcdf(outDir + '/' + outFile)


def process_GLODAP(inDir, outDir):
    """
    Unzip the gzipped data, trim folder, and process remaining data for
    MPAS-Analysis
    """
    # Authors
    # -------
    # Riley X. Brady

    # Check if already unzipped
    if os.path.isdir(inDir + '/GLODAPv2.2016b_MappedClimatologies'):
        print('GLODAP observations already unzipped.')
    # Unzip otherwise
    else:
        print('Unzipping GLODAP data...')
        tar = tarfile.open(inDir + '/GLODAPv2.2016b_MappedClimatologies.tar.gz')
        tar.extractall(inDir)
        tar.close()
    
    # Loop through files and delete if they aren't DIC, Alk
    print('Deleting unneeded observations...')
    keep_vars = ['TCO2', 'TAlk', 'pHtsinsitutp']
    for filename in os.listdir(inDir + '/GLODAPv2.2016b_MappedClimatologies'):
        if not any(s in filename for s in keep_vars):
            os.remove(inDir + '/GLODAPv2.2016b_MappedClimatologies/' + 
                      filename)
    
    # Edit to align with MPAS-Analysis standards
    updated_name = {'TCO2': 'DIC', 'TAlk': 'ALK', 'PI_TCO2': 'PI_DIC',
                    'pHtsinsitutp': 'pH_3D'}
    for v in ['TCO2', 'TAlk', 'PI_TCO2', 'pHtsinsitutp']:
        if not os.path.isfile(outDir + '/' + updated_name[v] + 
                              '_1.0x1.0degree.nc'):
            print("Processing and saving " + v + "...")
            filename = (inDir + '/GLODAPv2.2016b_MappedClimatologies/' + 
                        'GLODAPv2.2016b.' + v + '.nc')
            ds = xr.open_dataset(filename)
            ds = ds[v].isel(depth_surface=0)
            temp_vals = ds.values
            # Repeat the annual data into 12 months to fake the system
            temp_vals = np.repeat(temp_vals[np.newaxis, :, :], 12, axis=0)
            ds = xr.DataArray(temp_vals, dims=['Time', 'lat', 'lon'],
                             coords=[range(0, 12), ds.lat, ds.lon])
            ds.coords['month'] = ('Time', range(1, 13))
            ds.coords['year'] = ('Time', np.ones(12))
            ds.name = updated_name[v]
            outFile = (ds.name + '_1.0x1.0degree.nc')
            ds.to_dataset().to_netcdf(outDir + '/' + outFile)


def process_landschuetzer(inDir, outDir):
    """
    Edit Landschuetzerv2016 raw netCDF output to conform with MPAS-Analysis
    infrastructure
    """
    # Authors
    # -------
    # Riley X. Brady
    updated_name = {'fgco2': 'CO2_gas_flux', 'spco2': 'pCO2surface'}
    for v in ['fgco2', 'spco2']:
        if not os.path.isfile(outDir + '/' + updated_name[v] +
                              '_1.0x1.0degree.nc'):
            print("Processing and saving " + updated_name[v] + "...") 
            ds = xr.open_dataset(inDir + '/spco2_1982-2015_MPI_SOM-FFN_v2016.nc',
                             drop_variables='date')
            ds = ds[v + '_smoothed']
            ds = ds.rename({'time': 'Time'})
            ds.name = updated_name[v] 
            ds.coords['month'] = ds['Time.month']
            ds.coords['year'] = ds['Time.year']
            outFile = (updated_name[v] + '_1.0x1.0degree.nc')
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
   
    if '~' in args.inDir:
        raise ValueError("""
            Please avoid using '~' for your input directory.
            The tarfile package doesn't like it. Input the full path.
            """)

    # + + + SeaWIFS data + + +
    seawifs_files = ['S19972442010273.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19972742010304.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19973052010334.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19973352010365.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19980012010031.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19980322010059.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19980602010090.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19980912010120.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19981212010151.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19981522010181.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19981822010212.L3m_MC_CHL_chlor_a_9km.nc',
                     'S19982132010243.L3m_MC_CHL_chlor_a_9km.nc']
    urlBase = 'https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/'
    download_files(seawifs_files, urlBase=urlBase, outDir=args.inDir)
    process_SeaWIFS(args.inDir, args.outDir)

    # + + + GLODAPv2 data + + +
    urlBase = ('https://www.nodc.noaa.gov/archive/arc0107/0162565/2.2/data/' +
               '0-data/mapped/')
    download_files(['GLODAPv2.2016b_MappedClimatologies.tar.gz'],
                    urlBase=urlBase, outDir=args.inDir)
    process_GLODAP(args.inDir, args.outDir)

    # + + + Landschuetzer Carbon Flux + + +
    urlBase = ('https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/' +
               'data/0-data/')
    download_files(['spco2_1982-2015_MPI_SOM-FFN_v2016.nc'],
                   urlBase=urlBase, outDir=args.inDir)
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
     
