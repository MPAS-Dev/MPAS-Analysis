MPAS-Seaice
-----------

The Model for Prediction Across Scales Sea Ice (MPAS-Seaice)
is designed for the simulations of sea ice on unstructured grids supported by
the MPAS framework.  The model has not yet been publicly released and does not
have public documentation.

Setting up Standalone MPAS Sea Ice Runs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In order to support all sea=ice analysis tasks from MPAS-Analysis, certain
"analysis members", Fortran modules that perform analysis during the
simulation, need to be enabled.

The following is a list of suggested values for namelist options, typically
found in ``namelist.seaice`` or ``mpassi_in`` (or ``mpas-cice_in`` in
older E3SM runs)::

     config_AM_timeSeriesStatsMonthly_enable = .true.

Additionally, the duration of the run should be set to at least two years and
typically longer before most analysis is useful::

     config_run_duration = '0002-00-00_00:00:00'

Several streams must be defined in the streams file, typically
``streams.seaice`` or ``streams.cice`` in older E3SM runs, (even if they will
not be written out -- ``output_interval="none"``)::

  <stream name="timeSeriesStatsMonthlyRestart"
          type="input;output"
          io_type="pnetcdf"
          filename_template="mpasseaice.rst.am.timeSeriesStatsMonthly.$Y-$M-$D_$S.nc"
          filename_interval="output_interval"
          clobber_mode="truncate"
          packages="timeSeriesStatsMonthlyAMPKG"
          input_interval="initial_only"
          output_interval="stream:restart:output_interval">
  </stream>

  <stream name="timeSeriesStatsMonthlyOutput"
          type="output"
          io_type="pnetcdf"
          filename_template="mpasseaice.hist.am.timeSeriesStatsMonthly.$Y-$M-$D.nc"
          filename_interval="00-01-00_00:00:00"
          output_interval="00-01-00_00:00:00"
          clobber_mode="truncate"
          packages="timeSeriesStatsMonthlyAMPKG">

          <var name="icePresent"/>
          <var name="iceAreaCell"/>
          <var name="iceVolumeCell"/>
          <var name="xtime"/>
  </stream>

The ``filename_tempalate`` can be modified as desired (these are the defalult
values from E3SM).  For the ``timeSeriesStatsMonthlyOutput`` stream, both the
filename_interval and the output_interval must currently be monthly
(``"0000-01-00_00:00:00"``).

Additional fields can be included in the ``timeSeriesStatsMonthlyOutput``
streams.  These are the minimum that allow the analysis to run successfully.
