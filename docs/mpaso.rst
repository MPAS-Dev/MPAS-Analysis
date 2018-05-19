MPAS Ocean
----------

The Model for Prediction Across Scales Ocean (MPAS-O) is designed for the
simulation of the ocean system from time scales of months to millenia and
spatial scales from sub 1 km to global circulations.

MPAS-O has demonstrated the ability to accurately reproduce mesoscale ocean
activity with a local mesh refinement strategy.

In addition to faciliating the study of multiscale phenomena within the ocean
system, MPAS-O is intended for the study of anthropogenic climate change as
the ocean component of climate system models.


Full documentaiton is available at:
https://mpas-dev.github.io/ocean/ocean.html

Setting up Standalone MPAS-O Runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to support all ocean analysis tasks from MPAS-Analysis, certain
"analysis members", Fortran modules that perform analysis during the
simulation, need to be enabled.

The following is a list of suggested values for namelist options, typically
found in ``namelist.ocean`` or ``mpaso_in`` (or ``mpas-o_in`` in older E3SM runs)::

   config_AM_surfaceAreaWeightedAverages_enable = .true.
   config_AM_surfaceAreaWeightedAverages_compute_interval = '0000-00-00_01:00:00'
   config_AM_layerVolumeWeightedAverage_enable = .true.
   config_AM_layerVolumeWeightedAverage_compute_interval = '0000-00-00_01:00:00'
   config_AM_meridionalHeatTransport_enable = .true.
   config_AM_meridionalHeatTransport_compute_interval = '0000-00-00_01:00:00'
   config_AM_mixedLayerDepths_enable = .true.
   config_AM_timeSeriesStatsMonthly_enable = .true.

Additionally, the duration of the run should be set to at least two years and
typically longer before most analysis is useful::

   config_run_duration = '0002-00-00_00:00:00'

Several streams must be defined in the streams file, typically
``streams.ocean``, (even if they will not be written out --
``output_interval="none"``)::

  <stream name="timeSeriesStatsMonthlyRestart"
          type="input;output"
          filename_template="restarts/restart.AM.timeSeriesStatsMonthly.$Y-$M-$D_$h.$m.$s.nc"
          filename_interval="output_interval"
          reference_time="0001-01-01_00:00:00"
          clobber_mode="truncate"
          packages="timeSeriesStatsMonthlyAMPKG"
          input_interval="initial_only"
          output_interval="stream:restart:output_interval" >
  </stream>

  <stream name="timeSeriesStatsMonthlyOutput"
          type="output"
          filename_template="analysis_members/timeSeriesStatsMonthly.$Y-$M.nc"
          filename_interval="0000-01-00_00:00:00"
          reference_time="0001-01-01_00:00:00"
          clobber_mode="truncate"
          packages="timeSeriesStatsMonthlyAMPKG"
          output_interval="00-01-00_00:00:00" >

          <var_array name="activeTracers"/>
          <var name="normalVelocity"/>
          <var name="vertVelocityTop"/>
          <var_array name="avgValueWithinOceanRegion"/>
          <var_array name="avgValueWithinOceanLayerRegion"/>
          <var name="dThreshMLD"/>
          <var name="meridionalHeatTransportLatZ"/>
          <var name="meridionalHeatTransportLat"/>
          <var name="binBoundaryMerHeatTrans"/>
          <var name="xtime"/>
  </stream>

  <stream name="layerVolumeWeightedAverageOutput"
          type="output"
          io_type="pnetcdf"
          filename_template="mpaso.hist.am.layerVolumeWeightedAverage.$Y-$M-$D.nc"
          filename_interval="00-01-00_00:00:00"
          output_interval="none"
          clobber_mode="truncate"
          packages="layerVolumeWeightedAverageAMPKG">

      <var name="xtime"/>
      <var name="daysSinceStartOfSim"/>
      <var_array name="minValueWithinOceanLayerRegion"/>
      <var_array name="maxValueWithinOceanLayerRegion"/>
      <var_array name="avgValueWithinOceanLayerRegion"/>
      <var_array name="minValueWithinOceanVolumeRegion"/>
      <var_array name="maxValueWithinOceanVolumeRegion"/>
      <var_array name="avgValueWithinOceanVolumeRegion"/>
  </stream>

  <stream name="meridionalHeatTransportOutput"
          type="output"
          io_type="pnetcdf"
          filename_template="mpaso.hist.am.meridionalHeatTransport.$Y-$M-$D.nc"
          filename_interval="00-01-00_00:00:00"
          output_interval="none"
          clobber_mode="truncate"
          packages="meridionalHeatTransportAMPKG">

      <var name="xtime"/>
      <var name="daysSinceStartOfSim"/>
      <var name="binBoundaryMerHeatTrans"/>
      <var name="meridionalHeatTransportLatZ"/>
      <var name="meridionalHeatTransportLat"/>
      <var name="refZMid"/>
      <var name="refBottomDepth"/>
  </stream>

  <stream name="surfaceAreaWeightedAveragesOutput"
          type="output"
          io_type="netcdf"
          filename_template="mpaso.hist.am.surfaceAreaWeightedAverages.$Y-$M-$D.nc"
          filename_interval="00-01-00_00:00:00"
          output_interval="none"
          clobber_mode="truncate"
          packages="surfaceAreaWeightedAveragesAMPKG">

      <var name="xtime"/>
      <var name="daysSinceStartOfSim"/>
      <var_array name="minValueWithinOceanRegion"/>
      <var_array name="maxValueWithinOceanRegion"/>
      <var_array name="avgValueWithinOceanRegion"/>
  </stream>

  <stream name="mixedLayerDepthsOutput"
          type="output"
          io_type="pnetcdf"
          filename_template="mpaso.hist.am.mixedLayerDepths.$Y-$M-$D.nc"
          filename_interval="00-01-00_00:00:00"
          output_interval="none"
          clobber_mode="truncate"
          packages="mixedLayerDepthsAMPKG">

      <var name="xtime"/>
      <var name="daysSinceStartOfSim"/>
      <stream name="mesh"/>
      <var name="tThreshMLD"/>
      <var name="dThreshMLD"/>
      <var name="tGradMLD"/>
      <var name="dGradMLD"/>
  </stream>

The ``filename_tempalate`` can be modified as desired (in most cases, these are
the defalult values from E3SM).  For the ``timeSeriesStatsMonthlyOutput``
stream, both the filename_interval and the output_interval must currently be
monthly (``"0000-01-00_00:00:00"``).

Additional fields can be included in the ``timeSeriesStatsMonthlyOutput``
streams.  These are the minimum that allow the analysis to run successfully.
