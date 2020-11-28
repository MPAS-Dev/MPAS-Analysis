.. _task_timeSeriesTransport:

timeSeriesTransport
======================

An analysis task for plotting time series of of temperature, salinity,
potential density and other fields of interest averaged over ocean regions.

Component and Tags::

  component: ocean
  tags: timeSeries, transport

Configuration Options
---------------------

The following configuration options are available for this task::

    [timeSeriesTransport]
    ## options related to plotting time series of transport through transects

    # list of ocean transects from geometric_features to plot or ['all'] for all
    # available transects.
    transectsToPlot = ['Drake Passage', 'Tasmania-Ant', 'Africa-Ant', 'Antilles Inflow',
                       'Mona Passage', 'Windward Passage', 'Florida-Cuba', 'Florida-Bahamas',
                       'Indonesian Throughflow', 'Agulhas', 'Mozambique Channel', 'Bering Strait',
                       'Lancaster Sound', 'Fram Strait', 'Nares Strait']

    # Number of months over which to compute moving average
    movingAverageMonths = 1

    # An optional first year for the tick marks on the x axis. Leave commented out
    # to start at the beginning of the time series.

    # firstYearXTicks = 1

    # An optional number of years between tick marks on the x axis.  Leave
    # commented out to determine the distance between ticks automatically.

    # yearStrideXTicks = 1

Transect Names
--------------

The ``transectsToPlot`` can be set to ``['all']`` (the default) to plot all of the
transects defined in ``transportTransects20200621.geojson``.  These are::

  ["Africa-Ant",  "Agulhas", "Antarctic Peninsula", "Antilles Inflow",
   "Baja CA blockage", "Baltic Sea Deepen", "Barents Sea Opening",
   "Bering Strait", "Davis Strait", "Drake Passage", "English Channel Deepen",
   "Florida-Bahamas", "Florida-Cuba", "Fram Strait", "Indonesian Throughflow",
   "Ireland North Channel Deepen", "Japan Hokkaido blockage",
   "Japan La Perouse Strait Deepen",  "Japan Tsugaru Strait Deepen",
   "Japan blockage", "Lancaster Sound", "Mona Passage", "Mozambique Channel",
   "Nares Strait", "Nares Strait Deepen", "Persian Gulf Deepen",
   "Red Sea Deepen", "Sakhalin blockage", "Strait of Gibralter Deepen 1",
    "Strait of Gibralter Deepen 2", "Tasmania-Ant", "White Sea",
    "Windward Passage"]

Many of these are likely not of interest in most simulations, so a subset of
the most relevant transects has been chosen in the default configuration.

Other Options
-------------

* :ref:`config_moving_average`
* :ref:`config_time_axis_ticks`

Example Result
--------------

.. image:: examples/drake_passage_transport.png
   :width: 500 px
   :align: center

.. _`antarctic_ocean_regions`: https://github.com/MPAS-Dev/geometric_features/tree/master/feature_creation_scripts/antarctic_ocean_regions
