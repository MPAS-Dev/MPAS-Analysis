'''
Mappings of namelist options, stream names and variable names from various
MPAS-SeaIce versions to those used by mpas_analysis

Authors
-------
Xylar Asay-Davis

Last Modified
-------------
03/29/2017
'''

seaIceNamelistMap = {}

seaIceStreamMap = {'timeSeriesStats': ['timeSeriesStatsMonthlyOutput']}


seaIceVariableMap = {}
seaIceVariableMap['Time'] = \
    [['xtime_startMonthly', 'xtime_endMonthly'],
     ['xtime_start', 'xtime_end'],
     'timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1',
     'time_avg_daysSinceStartOfSim']

seaIceVariableMap['iceAreaCell'] = \
    ['timeMonthly_avg_iceAreaCell',
     'timeSeriesStatsMonthly_avg_iceAreaCell_1']

seaIceVariableMap['iceVolumeCell'] = \
    ['timeMonthly_avg_iceVolumeCell',
     'timeSeriesStatsMonthly_avg_iceVolumeCell_1']
