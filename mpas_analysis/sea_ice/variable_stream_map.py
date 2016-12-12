# mappings of stream names from various MPAS-SI versions to those in
# mpas_analysis
seaIceStreamMap = {'timeSeriesStats': ['timeSeriesStatsMonthlyOutput']}


seaIceVariableMap = {}
seaIceVariableMap['Time'] = \
    [['xtime_start', 'xtime_end'],
     'timeSeriesStatsMonthly_avg_daysSinceStartOfSim_1',
     'time_avg_daysSinceStartOfSim']

seaIceVariableMap['iceAreaCell'] = \
    ['timeSeriesStatsMonthly_avg_iceAreaCell_1']

seaIceVariableMap['iceVolumeCell'] = \
    ['timeSeriesStatsMonthly_avg_iceVolumeCell_1']
