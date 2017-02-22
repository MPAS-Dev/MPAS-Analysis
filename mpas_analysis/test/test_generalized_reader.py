"""
Unit test infrastructure for the generalized_reader.

Xylar Asay-Davis
02/16/2017
"""

import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.generalized_reader.generalized_reader \
    import open_multifile_dataset


@pytest.mark.usefixtures("loaddatadir")
class TestGeneralizedReader(TestCase):

    def test_variableMap(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        simulationStartTime = '0001-01-01'
        variableMap = {
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

        variableList = ['avgSurfaceTemperature', 'avgLayerTemperature',
                        'refBottomDepth', 'daysSinceStartOfSim']

        for calendar in ['gregorian', 'gregorian_noleap']:
            # preprocess_mpas will use variableMap to map the variable names
            # from their values in the file to the desired values in
            # variableList
            ds = open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                simulationStartTime=simulationStartTime,
                timeVariableName='Time',
                variableList=variableList,
                variableMap=variableMap)

            # make sure the remapping happened as expected
            self.assertEqual(sorted(ds.data_vars.keys()), sorted(variableList))

    def test_open_dataset_fn(self):
        fileName = str(self.datadir.join('example_jan.nc'))
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        for calendar in ['gregorian', 'gregorian_noleap']:
            ds = open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                timeVariableName=timestr,
                variableList=variableList)
            self.assertEqual(ds.data_vars.keys(), variableList)

    def test_start_end(self):
        fileName = str(self.datadir.join('example_jan_feb.nc'))
        timestr = ['xtime_start', 'xtime_end']
        variableList = \
            ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

        for calendar in ['gregorian', 'gregorian_noleap']:
            # all dates
            ds = open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                timeVariableName=timestr,
                variableList=variableList,
                startDate='0001-01-01',
                endDate='9999-12-31')
            self.assertEqual(len(ds.Time), 2)

            # just the first date
            ds = open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                timeVariableName=timestr,
                variableList=variableList,
                startDate='0005-01-01',
                endDate='0005-02-01')
            self.assertEqual(len(ds.Time), 1)

            # just the second date
            ds = open_multifile_dataset(
                fileNames=fileName,
                calendar=calendar,
                timeVariableName=timestr,
                variableList=variableList,
                startDate='0005-02-01',
                endDate='0005-03-01')
            self.assertEqual(len(ds.Time), 1)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
