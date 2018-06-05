# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#
import numpy
import xarray
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import pyproj
import sys

from mpas_analysis.shared.interpolation import Remapper
from mpas_analysis.shared.grid import ProjectionGridDescriptor
from mpas_analysis.shared.mpas_xarray.mpas_xarray import subset_variables
from mpas_analysis.shared.climatology \
    import get_Antarctic_stereographic_comparison_descriptor
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser

inFileName = '/media/xylar/extra_data/data_overflow/observations/Antarctica/' \
             'Rignot_et_al._2013/Ant_MeltingRate.nc'

config = MpasAnalysisConfigParser()
config.read('config.default')

ds = xarray.open_dataset(inFileName)
ds = subset_variables(ds, ['melt_actual', 'xaxis', 'yaxis'])
lx = numpy.abs(1e-3*(ds.xaxis.values[-1]-ds.xaxis.values[0]))
ly = numpy.abs(1e-3*(ds.yaxis.values[-1]-ds.yaxis.values[0]))

maskedMeltRate = numpy.ma.masked_array(ds.melt_actual,
                                       mask=(ds.melt_actual.values == 0.))

ds['meltRate'] = xarray.DataArray(maskedMeltRate, dims=ds.melt_actual.dims,
                                  coords=ds.melt_actual.coords,
                                  attrs=ds.melt_actual.attrs)

ds = ds.drop('melt_actual')

inGridName = '{}x{}km_1.0km_Antarctic_stereo'.format(lx, ly)

projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                         '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

inDescriptor = ProjectionGridDescriptor(projection)

inDescriptor.read(inFileName,  xVarName='xaxis', yVarName='yaxis',
                  meshName=inGridName)

outDescriptor = get_Antarctic_stereographic_comparison_descriptor(config)
outGridName = outDescriptor.meshName

outFileName = 'Rignot_2013_melt_rates_{}.nc'.format(outGridName)

mappingFileName = 'map_{}_to_{}.nc'.format(inGridName, outGridName)

remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

remapper.build_mapping_file(method='bilinear')

remappedDataset = remapper.remap(ds, renormalizationThreshold=0.01)

remappedDataset.attrs['history'] = ' '.join(sys.argv)
remappedDataset.to_netcdf(outFileName)

norm = colors.SymLogNorm(linthresh=1, linscale=1, vmin=-100.0, vmax=100.0)

plt.figure()
plt.imshow(maskedMeltRate, origin='upper', norm=norm)
plt.colorbar()
plt.figure()
plt.imshow(remappedDataset.meltRate.values, origin='lower', norm=norm)
plt.colorbar()

plt.show()
