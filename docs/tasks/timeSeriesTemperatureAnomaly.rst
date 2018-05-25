.. _task_timeSeriesTemperatureAnomaly:

timeSeriesTemperatureAnomaly
============================

An analysis task for plotting a Hovmoller plot (time and depth axes) of the
anomaly in ocean potential temperature from a reference year (usully the first
year of the simulation).

Component and Tags::

  component: ocean
  tags: timeSeries, temperature, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [hovmollerTemperatureAnomaly]
  ## options related to plotting time series of potential temperature vs. depth

  # list of regions to plot from the region list in [regions] below
  regions = ['global']

  # Number of points over which to compute moving average(e.g., for monthly
  # output, movingAveragePoints=12 corresponds to a 12-month moving average
  # window)
  movingAveragePoints = 12

  # colormap
  colormapName = balance
  # color indices into colormapName for filled contours
  colormapIndices = [0, 28, 57, 85, 113, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevels = [-1, -0.5, -0.2, -0.05, 0, 0.05, 0.2, 0.5, 1]
  # contour line levels
  contourLevels = np.arange(-1.0, 1.26, 0.25)

  # An optional first year for the tick marks on the x axis. Leave commented out
  # to start at the beginning of the time series.

  # firstYearXTicks = 1

  # An optional number of years between tick marks on the x axis.  Leave
  # commented out to determine the distance between ticks automatically.

  # yearStrideXTicks = 1

For more details, see:
 * :ref:`config_regions`
 * :ref:`config_colormaps`
 * :ref:`config_moving_average`
 * :ref:`config_time_axis_ticks`

Example Result
--------------

.. image:: examples/hovmoller_temp.png
   :width: 500 px
   :align: center
