# mappings of stream names from various MPAS-O versions to those in
# mpas_analysis
oceanStreamMap = {'timeSeriesStats': ['timeSeriesStatsOutput',
                                      'timeSeriesStatsMonthly',
                                      'timeSeriesStatsMonthlyOutput']}


oceanVariableMap = {}
oceanVariableMap['Time'] = [['xtime_startMonthly', 'xtime_endMonthly'],
                            ['xtime_start', 'xtime_end'],
                            'time_avg_daysSinceStartOfSim']

# SST
oceanVariableMap['avgSurfaceTemperature'] = \
    ['time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature',
     'time_avg_avgValueWithinOceanRegion_avgSurfaceTemperature_1',
     'timeMonthly_avg_avgValueWithinOceanRegion_avgSurfaceTemperature']

# OHC
oceanVariableMap['avgLayerTemperature'] = \
    ['time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature',
     'time_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature_1',
     'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerTemperature']
oceanVariableMap['sumLayerMaskValue'] = \
    ['time_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue',
     'time_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue_1',
     'timeMonthly_avg_avgValueWithinOceanLayerRegion_sumLayerMaskValue']
oceanVariableMap['avgLayerArea'] = \
    ['time_avg_avgValueWithinOceanLayerRegion_avgLayerArea',
     'time_avg_avgValueWithinOceanLayerRegion_avgLayerArea_1',
     'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerArea']
oceanVariableMap['avgLayerThickness'] = \
    ['time_avg_avgValueWithinOceanLayerRegion_avgLayerThickness',
     'time_avg_avgValueWithinOceanLayerRegion_avgLayerThickness_1',
     'timeMonthly_avg_avgValueWithinOceanLayerRegion_avgLayerThickness']

# MOC
oceanVariableMap['avgNormalVelocity'] = \
    ['time_avg_normalVelocity',
     'time_avg_normalVelocity_1',
     'timeMonthly_avg_normalVelocity']
oceanVariableMap['avgVertVelocityTop'] = \
    ['time_avg_vertVelocityTop',
     'time_avg_vertVelocityTop_1',
     'timeMonthly_avg_vertVelocityTop']

# model vs. obs.
oceanVariableMap['mld'] = \
    ['time_avg_dThreshMLD',
     'time_avg_dThreshMLD_1',
     'timeMonthly_avg_dThreshMLD']

oceanVariableMap['sst'] = \
    ['time_avg_activeTracers_temperature',
     'time_avg_activeTracers_temperature_1',
     'timeMonthly_avg_activeTracers_temperature']

oceanVariableMap['sss'] = \
    ['time_avg_activeTracers_salinity',
     'time_avg_activeTracers_salinity_1',
     'timeMonthly_avg_activeTracers_salinity']
