"""
Unit test infrastructure for horizontal interpolation.

Xylar Asay-Davis
02/25/2017
"""

import pytest
import shutil
import os
import tempfile
import numpy

from mpas_analysis.shared.interpolation import interpolate
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


@pytest.mark.usefixtures("loaddatadir")
class TestInterp(TestCase):

    def setUp(self):
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.test_dir)

    def test_destination_from_grid_file(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from a file containing 'lat' and 'lon' coords

        Xylar Asay-Davis
        02/25/2017
        """

        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        latLonGridFileName = str(self.datadir.join('SST_annual_1870-1900.nc'))
        timeSeriesFileName = str(self.datadir.join('timeSeries.0002-01-01.nc'))

        suffix = 'destination_from_grid_file'
        weightFileName = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outFileName = '{}/remapped_{}.nc'.format(self.test_dir, suffix)

        interpolate.build_remap_weights(sourceFileName=mpasMeshFileName,
                                        outWeightFileName=weightFileName,
                                        destintionFileName=latLonGridFileName,
                                        destintionLatVarName='lat',
                                        destintionLonVarName='lon',
                                        sourceFileType='mpas',
                                        method='bilinear')

        assert os.path.exists(weightFileName)

        interpolate.remap(inFileName=timeSeriesFileName,
                          outFileName=outFileName,
                          inWeightFileName=weightFileName,
                          sourceFileType='mpas')

        assert os.path.exists(outFileName)

        # TODO: check the results against a reference result

    def test_destination_from_numpy_lat_lon(self):
        """
        test horizontal interpolation from an MPAS mesh to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.

        Xylar Asay-Davis
        02/25/2017
        """

        configPath = str(self.datadir.join('config.analysis'))
        config = MpasAnalysisConfigParser()
        config.read(configPath)

        lat = numpy.array(config.getExpression('interpolate', 'lat',
                                               usenumpyfunc=True))
        lon = numpy.array(config.getExpression('interpolate', 'lon',
                                               usenumpyfunc=True))

        mpasMeshFileName = str(self.datadir.join('mpasMesh.nc'))
        timeSeriesFileName = str(self.datadir.join('timeSeries.0002-01-01.nc'))

        suffix = 'destination_from_config_options'
        weightFileName = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outFileName = '{}/remapped_{}.nc'.format(self.test_dir, suffix)

        interpolate.build_remap_weights(sourceFileName=mpasMeshFileName,
                                        outWeightFileName=weightFileName,
                                        sourceFileType='mpas',
                                        method='bilinear',
                                        destinationLat=lat,
                                        destinationLon=lon)

        assert os.path.exists(weightFileName)

        interpolate.remap(inFileName=timeSeriesFileName,
                          outFileName=outFileName,
                          inWeightFileName=weightFileName,
                          sourceFileType='mpas')

        assert os.path.exists(outFileName)

        # TODO: check the results against a reference result

    def test_source_lat_lon(self):
        """
        test horizontal interpolation from a lat/lon grid to a destination
        lat/lon grid determined from config options 'lat' and 'lon'.

        Xylar Asay-Davis
        02/25/2017
        """

        lat = numpy.linspace(-90., 90., 361)
        lon = numpy.linspace(-180., 180., 721)

        sourceFileName = str(self.datadir.join('SST_annual_1870-1900.nc'))

        suffix = 'source_lat_lon'
        weightFileName = '{}/weights_{}.nc'.format(self.test_dir, suffix)
        outFileName = '{}/remapped_{}.nc'.format(self.test_dir, suffix)

        interpolate.build_remap_weights(sourceFileName=sourceFileName,
                                        outWeightFileName=weightFileName,
                                        sourceFileType='latlon',
                                        method='bilinear',
                                        destinationLat=lat,
                                        destinationLon=lon)

        assert os.path.exists(weightFileName)

        interpolate.remap(inFileName=sourceFileName,
                          outFileName=outFileName,
                          inWeightFileName=weightFileName,
                          sourceFileType='mpas')

        assert os.path.exists(outFileName)

        # TODO: check the results against a reference result

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
