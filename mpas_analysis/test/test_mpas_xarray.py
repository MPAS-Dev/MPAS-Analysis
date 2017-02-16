"""
Unit test infrastructure for mpas_xarray.

Xylar Asay-Davis, Phillip J. Wolfram
02/15/2017
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.mpas_xarray import mpas_xarray
import pandas as pd


@pytest.mark.usefixtures("loaddatadir")
class TestMpasXarray(TestCase):

    def test_subset_variables(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        # first, test loading the whole data set and then calling
        # subset_variables explicitly
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                timeVariableName=timestr,
                                                yearOffset=1850)
        ds = mpas_xarray.subset_variables(ds, variableList)
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        self.assertEqual(pd.Timestamp(ds.Time.values[0]),
                         pd.Timestamp('1855-01-16 12:22:30'))

        # next, test the same with the onlyvars argument
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                timeVariableName=timestr,
                                                variableList=variableList,
                                                yearOffset=1850)
        self.assertEqual(ds.data_vars.keys(), variableList)

        with self.assertRaisesRegexp(ValueError,
                                     'Empty dataset is returned.'):
            missingvars = ['foo', 'bar']
            ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                    timeVariableName=timestr,
                                                    variableList=missingvars,
                                                    yearOffset=1850)

    def test_iselvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        iselvals = {'nVertLevels': slice(0, 3)}
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                timeVariableName=timestr,
                                                variableList=variableList,
                                                iselValues=iselvals,
                                                yearOffset=1850)

        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        self.assertEqual(ds[variableList[0]].shape, (1, 7, 3))
        self.assertEqual(ds['refBottomDepth'].shape, (3,))
        self.assertApproxEqual(ds['refBottomDepth'][-1],
                               4.882000207901)
        date = pd.Timestamp(ds.Time.values[0])
        # round to nearest second
        date = pd.Timestamp(long(round(date.value, -9)))
        self.assertEqual(date, pd.Timestamp('1855-01-13 12:24:14'))

    def test_no_units(self):
        fileName = str(self.datadir.join('example_no_units_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                timeVariableName=timestr,
                                                variableList=variableList,
                                                yearOffset=1850)
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        date = pd.Timestamp(ds.Time.values[0])
        # round to nearest second
        date = pd.Timestamp(long(round(date.value, -9)))
        self.assertEqual(date, pd.Timestamp('1855-01-13 12:24:14'))

    def test_bad_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        selvals = {'refBottomDepth': 8.77999997138977}
        with self.assertRaisesRegexp(AssertionError,
                                     'not a dimension in the dataset that '
                                     'can be used for selection'):
            mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                               timeVariableName=timestr,
                                               variableList=variableList,
                                               selValues=selvals,
                                               yearOffset=1850)

    def test_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        dsRef = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                   timeVariableName=timestr,
                                                   variableList=variableList,
                                                   selValues=None,
                                                   yearOffset=1850)

        for vertIndex in range(0, 11):
            selvals = {'nVertLevels': vertIndex}
            ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                    timeVariableName=timestr,
                                                    variableList=variableList,
                                                    selValues=selvals,
                                                    yearOffset=1850)

            self.assertEqual(ds.data_vars.keys(), variableList)
            self.assertEqual(ds[variableList[0]].shape, (1, 7))
            self.assertEqual(ds['refBottomDepth'],
                             dsRef['refBottomDepth'][vertIndex])

    def test_remove_repeated_time_index(self):
        fileName = str(self.datadir.join('example_jan*.nc'))
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        # repeat time indices are removed in openMultifileDataSet
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                timeVariableName=timestr,
                                                variableList=variableList,
                                                yearOffset=1850)

        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        # There would be 3 time indices if repeat indices had not been removed.
        # Make sure there are 2.
        self.assertEqual(len(ds.Time.values), 2)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
