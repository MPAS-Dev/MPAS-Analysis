"""
Unit test infrastructure for mpas_xarray.

Xylar Asay-Davis, Phillip J. Wolfram
02/15/2017
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.mpas_xarray import mpas_xarray
from mpas_analysis.shared.timekeeping.utility import days_to_datetime, \
    string_to_datetime


@pytest.mark.usefixtures("loaddatadir")
class TestMpasXarray(TestCase):

    def test_subset_variables(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        calendar = 'gregorian_noleap'
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        # first, test loading the whole data set and then calling
        # subset_variables explicitly
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                calendar=calendar,
                                                timeVariableName=timestr)
        ds = mpas_xarray.subset_variables(ds, variableList)
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        self.assertEqual(days_to_datetime(days=ds.Time.values,
                                          referenceDate='0001-01-01',
                                          calendar=calendar),
                         string_to_datetime('0005-01-16 12:22:30'))
        # next, test the same with the onlyvars argument
        ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                calendar=calendar,
                                                timeVariableName=timestr,
                                                variableList=variableList)
        self.assertEqual(ds.data_vars.keys(), variableList)

        with self.assertRaisesRegexp(ValueError,
                                     'Empty dataset is returned.'):
            missingvars = ['foo', 'bar']
            ds = mpas_xarray.open_multifile_dataset(fileNames=fileName,
                                                    calendar=calendar,
                                                    timeVariableName=timestr,
                                                    variableList=missingvars)

    def test_iselvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        calendar = 'gregorian_noleap'
        simulationStartTime = '0001-01-01'
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        iselvals = {'nVertLevels': slice(0, 3)}
        ds = mpas_xarray.open_multifile_dataset(
            fileNames=fileName,
            calendar=calendar,
            simulationStartTime=simulationStartTime,
            timeVariableName=timestr,
            variableList=variableList,
            iselValues=iselvals)

        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        self.assertEqual(ds[variableList[0]].shape, (1, 7, 3))
        self.assertEqual(ds['refBottomDepth'].shape, (3,))
        self.assertApproxEqual(ds['refBottomDepth'][-1],
                               4.882000207901)

        self.assertEqual(days_to_datetime(days=ds.Time.values[0],
                                          referenceDate='0001-01-01',
                                          calendar=calendar),
                         string_to_datetime('0005-01-14 12:24:14'))

    def test_no_units(self):
        fileName = str(self.datadir.join('example_no_units_jan.nc'))
        calendar = 'gregorian_noleap'
        simulationStartTime = '0001-01-01'
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        ds = mpas_xarray.open_multifile_dataset(
            fileNames=fileName,
            calendar=calendar,
            simulationStartTime=simulationStartTime,
            timeVariableName=timestr,
            variableList=variableList)
        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))

        self.assertEqual(days_to_datetime(days=ds.Time.values[0],
                                          referenceDate='0001-01-01',
                                          calendar=calendar),
                         string_to_datetime('0005-01-14 12:24:14'))

    def test_bad_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        calendar = 'gregorian_noleap'
        simulationStartTime = '0001-01-01'
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        selvals = {'refBottomDepth': 8.77999997138977}
        with self.assertRaisesRegexp(AssertionError,
                                     'not a dimension in the dataset that '
                                     'can be used for selection'):
            mpas_xarray.open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                simulationStartTime=simulationStartTime,
                timeVariableName=timestr,
                variableList=variableList,
                selValues=selvals)

    def test_selvals(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        calendar = 'gregorian_noleap'
        simulationStartTime = '0001-01-01'
        timestr = 'time_avg_daysSinceStartOfSim'
        variableList = \
            ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
             'refBottomDepth']

        dsRef = mpas_xarray.open_multifile_dataset(
            fileNames=fileName,
            calendar=calendar,
            simulationStartTime=simulationStartTime,
            timeVariableName=timestr,
            variableList=variableList,
            selValues=None)

        for vertIndex in range(0, 11):
            selvals = {'nVertLevels': vertIndex}
            ds = mpas_xarray.open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                simulationStartTime=simulationStartTime,
                timeVariableName=timestr,
                variableList=variableList,
                selValues=selvals)

            self.assertEqual(ds.data_vars.keys(), variableList)
            self.assertEqual(ds[variableList[0]].shape, (1, 7))
            self.assertEqual(ds['refBottomDepth'],
                             dsRef['refBottomDepth'][vertIndex])

    def test_remove_repeated_time_index(self):
        fileName = str(self.datadir.join('example_jan*.nc'))
        calendar = 'gregorian_noleap'
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        # repeat time indices are removed in openMultifileDataSet
        ds = mpas_xarray.open_multifile_dataset(
            fileNames=fileName,
            calendar=calendar,
            timeVariableName=timestr,
            variableList=variableList)

        self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))
        # There would be 3 time indices if repeat indices had not been removed.
        # Make sure there are 2.
        self.assertEqual(len(ds.Time.values), 2)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
