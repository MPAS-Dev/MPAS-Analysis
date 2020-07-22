.. _config_time_series:

Time Series
===========

The ``[timeSeries]`` section of a configuration file contains options used to
determine the start and end years of time series plots and the reference years
for anomalies::

  [timeSeries]
  ## options related to producing time series plots, often to compare against
  ## observations and previous runs

  # the year from which to compute anomalies if not the start year of the
  # simulation.  This might be useful if a long spin-up cycle is performed and
  # only the anomaly over a later span of years is of interest.
  # anomalyRefYear = 249

  # start and end years for timeseries analysis. Use endYear = end to indicate
  # that the full range of the data should be used.  If errorOnMissing = False,
  # the start and end year will be clipped to the valid range.  Otherwise, out
  # of bounds values will lead to an error.  In a "control" config file used in
  # a "main vs. control" analysis run, the range of years must be valid and
  # cannot include "end" because the original data may not be available.
  startYear = 1
  endYear = end

Start and End Year
------------------

A custom config file should specify a start and end year for time series.
If ``errorOnMissing = False`` in the ``input`` section and the start or end
year is beyond the range of the simulation, the range will be reduced to those
dates with available data and a warning message will be displayed.  If
``errorOnMissing = True``, out of range year will produce an error.


Anomaly Reference Year
----------------------

Anomalies between a moving average of a time series and the average over a
reference year are used in several analysis tasks.  By default, the reference
year is not specified in the configuration file and it taken to be the start of
the simulation (determined from the contents of a restart file).  Under certain
circumstances (e.g. repetition of forcing data for several cycles, as in
the `Common Ocean Reference Experiments, CORE`_), it may be desirable to
specify a different year to use for computing anomalies::

  anomalyRefYear = 249

.. _`Common Ocean Reference Experiments, CORE`: http://data1.gfdl.noaa.gov/nomads/forms/mom4/CORE.html
