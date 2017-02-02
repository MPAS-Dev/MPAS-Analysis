"""
Unit test infrastructure for mpas_xarray.

Xylar Asay-Davis, Phillip J. Wolfram
12/07/2016
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.mpas_xarray import mpas_xarray
import xarray as xr
import pandas as pd


@pytest.mark.usefixtures("loaddatadir")
class TestNamelist(TestCase):

    def test_subset_variables(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = ['xtime_start', 'xtime_end']
        varList = ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        # first, test loading the whole data set and then calling
        # subset_variables explicitly
        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             yearoffset=1850))
        ds = mpas_xarray.subset_variables(ds, varList)
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(varList))
        self.assertEqual(pd.Timestamp(ds.Time.values[0]),
                         pd.Timestamp('1855-01-16 12:22:30'))

        # next, test the same with the onlyvars argument
        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             onlyvars=varList,
                                                             yearoffset=1850))
        self.assertEqual(ds.data_vars.keys(), varList)

        with self.assertRaisesRegexp(AssertionError,
                                     'Empty dataset is returned.'):
            missingvars = ['foo', 'bar']
            ds = xr.open_mfdataset(
                fileName,
                preprocess=lambda x:
                    mpas_xarray.preprocess_mpas(x,
                                                timestr=timestr,
                                                onlyvars=missingvars,
                                                yearoffset=1850))

    def test_iselvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        varList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        iselvals = {'nVertLevels': slice(0, 3)}
        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             onlyvars=varList,
                                                             iselvals=iselvals,
                                                             yearoffset=1850))
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(varList))
        self.assertEqual(ds[varList[0]].shape, (1, 7, 3))
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
        varList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             onlyvars=varList,
                                                             yearoffset=1850))
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(varList))
        date = pd.Timestamp(ds.Time.values[0])
        # round to nearest second
        date = pd.Timestamp(long(round(date.value, -9)))
        self.assertEqual(date, pd.Timestamp('1855-01-13 12:24:14'))

    def test_bad_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        varList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        selvals = {'refBottomDepth': 8.77999997138977}
        with self.assertRaisesRegexp(AssertionError,
                                     'not a dimension in the dataset that '
                                     'can be used for selection'):
            ds = xr.open_mfdataset(
                fileName,
                preprocess=lambda x:
                    mpas_xarray.preprocess_mpas(x,
                                                timestr=timestr,
                                                onlyvars=varList,
                                                selvals=selvals,
                                                yearoffset=1850))

    def test_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = 'time_avg_daysSinceStartOfSim'
        varList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        dsRef = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             onlyvars=varList,
                                                             selvals=None,
                                                             yearoffset=1850))

        for vertIndex in range(0, 11):
            selvals = {'nVertLevels': vertIndex}
            ds = xr.open_mfdataset(
                fileName,
                preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                                 timestr=timestr,
                                                                 onlyvars=varList,
                                                                 selvals=selvals,
                                                                 yearoffset=1850))
            self.assertEqual(ds.data_vars.keys(), varList)
            self.assertEqual(ds[varList[0]].shape, (1, 7))
            self.assertEqual(ds['refBottomDepth'],
                             dsRef['refBottomDepth'][vertIndex])

    def test_remove_repeated_time_index(self):
        fileName = str(self.datadir.join('example_jan*.nc'))
        timestr = ['xtime_start', 'xtime_end']
        varList = ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(x,
                                                             timestr=timestr,
                                                             onlyvars=varList,
                                                             yearoffset=1850))

        self.assertEqual(sorted(ds.data_vars.keys()), sorted(varList))
        self.assertEqual(len(ds.Time.values), 3)

        ds = mpas_xarray.remove_repeated_time_index(ds)
        self.assertEqual(len(ds.Time.values), 2)

    def test_variable_map(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        varMap = {
            'avgSurfaceTemperature':
                ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature',
                 'other_string',
                 'yet_another_string'],
            'daysSinceStartOfSim':
                ['time_avg_daysSinceStartOfSim',
                 'xtime',
                 'something_else'],
            'avgLayerTemperature':
                ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
                 'test1',
                 'test2'],
            'Time': [['xtime_start', 'xtime_end'],
                     'time_avg_daysSinceStartOfSim']}

        varList = ['avgSurfaceTemperature', 'avgLayerTemperature',
                   'refBottomDepth', 'daysSinceStartOfSim']

        # preprocess_mpas will use varMap to map the variable names from their
        # values in the file to the desired values in varList
        ds = xr.open_mfdataset(
            fileName,
            preprocess=lambda x: mpas_xarray.preprocess_mpas(
                x,
                timestr='Time',
                onlyvars=varList,
                yearoffset=1850,
                varmap=varMap))

        # make sure the remapping happened as expected
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(varList))

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
