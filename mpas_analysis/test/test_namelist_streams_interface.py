"""
Unit test infrastructure for namelist and streams readers, adapted from
approach of xarray.

Phillip J. Wolfram, Xylar Asay-Davis
10/26/2016
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.io import NameList, StreamsFile


@pytest.mark.usefixtures("loaddatadir")
class TestNamelist(TestCase):
    def setup_namelist(self):
        nlpath = self.datadir.join('namelist.ocean')
        self.nl = NameList(bytes(nlpath))

    def setup_streams(self):
        sfpath = self.datadir.join('streams.ocean')
        self.sf = StreamsFile(bytes(sfpath))

    def test_open_files(self):
        self.setup_namelist()
        self.setup_streams()

    def test_read_namelist(self):
        self.setup_namelist()

        # check accessing generalized function techniques
        self.assertEqual(self.nl.config_dt, '00:10:00')
        self.assertEqual(self.nl['config_dt'], '00:10:00')

        # check cast accessors
        self.assertEqual(self.nl.getint('config_num_halos'), 3)
        self.assertApproxEqual(self.nl.getfloat('config_min_thickness'), 1.0)
        self.assertEqual(self.nl.getbool('config_do_restart'), False)

        # tests for use of ' and " for string selections
        self.assertEqual(self.nl.config_test_extra_equals1, 'a = b')
        self.assertEqual(self.nl.config_test_extra_equals2, 'a = b')

    def test_read_streamsfile(self):
        self.setup_streams()

        # check
        self.assertEqual(self.sf.read('output', 'type'), 'output')
        self.assertEqual(self.sf.read('restart', 'output_interval'),
                         '0100_00:00:00')

        files = self.sf.readpath('output')
        expectedFiles = []
        for date in ['0001-01-01', '0001-01-02', '0001-02-01', '0002-01-01']:
            expectedFiles.append('{}/output/output.{}_00.00.00.nc'
                                 .format(self.sf.streamsdir, date))
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('output',
                                 startDate='0001-01-03',
                                 endDate='0001-12-30',
                                 calendar='gregorian_noleap')
        expectedFiles = []
        for date in ['0001-01-02', '0001-02-01']:
            expectedFiles.append('{}/output/output.{}_00.00.00.nc'
                                 .format(self.sf.streamsdir, date))
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('output',
                                 startDate='0001-01-03',
                                 calendar='gregorian_noleap')
        expectedFiles = []
        for date in ['0001-01-02', '0001-02-01', '0002-01-01']:
            expectedFiles.append('{}/output/output.{}_00.00.00.nc'
                                 .format(self.sf.streamsdir, date))
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('output',
                                 endDate='0001-12-30',
                                 calendar='gregorian_noleap')
        expectedFiles = []
        for date in ['0001-01-01', '0001-01-02', '0001-02-01']:
            expectedFiles.append('{}/output/output.{}_00.00.00.nc'
                                 .format(self.sf.streamsdir, date))
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('restart',
                                 startDate='0001-01-01',
                                 endDate='0001-12-31',
                                 calendar='gregorian_noleap')
        expectedFiles = []
        for seconds in ['00010', '00020']:
            expectedFiles.append('{}/restarts/restart.0001-01-01_{}.nc'
                                 .format(self.sf.streamsdir, seconds))
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('mesh')
        expectedFiles = ['{}/mesh.nc'.format(self.sf.streamsdir)]
        self.assertEqual(files, expectedFiles)

        files = self.sf.readpath('mesh',
                                 startDate='0001-01-01',
                                 endDate='0001-12-31',
                                 calendar='gregorian_noleap')
        expectedFiles = ['{}/mesh.nc'.format(self.sf.streamsdir)]
        self.assertEqual(files, expectedFiles)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
