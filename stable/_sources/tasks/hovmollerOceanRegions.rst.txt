.. _task_oceanRegionalProfiles:

hovmollerOceanRegions
=====================

An analysis task for plotting depth profiles vs. time of temperature, salinity,
potential density, etc. averaged over regions.

Component and Tags::

  component: ocean
  tags: profiles, climatology, hovmoller

Configuration Options
---------------------

The following configuration options are available for this task:

.. code-block:: cfg

    [hovmollerOceanRegions]
    ## options related to plotting Hovmoller diagrams (depth vs. time plots) of
    ## regional means of 3D MPAS fields

    # the names of region groups to plot, each with its own section below
    regionGroups = ['Ocean Basins']


    [hovmollerOceanBasins]
    ## options related to plotting Hovmoller diagrams of ocean basins

    # a list of dictionaries for each field to plot.  The dictionary includes
    # prefix (used for file names, task names and sections) as well as the MPAS
    # name of the field, units for colorbars and a the name as it should appear
    # in figure titles and captions.
    fields =
        [{'prefix': 'potentialTemperature',
          'mpas': 'timeMonthly_avg_activeTracers_temperature',
          'units': r'$\degree$C',
          'titleName': 'Potential Temperature'},
         {'prefix': 'salinity',
          'mpas': 'timeMonthly_avg_activeTracers_salinity',
          'units': r'PSU',
          'titleName': 'Salinity'},
         {'prefix': 'potentialDensity',
          'mpas': 'timeMonthly_avg_potentialDensity',
          'units': r'kg m$^{-3}$',
          'titleName': 'Potential Density'}]

    # a list of region names from the region masks file to plot
    regionNames = ["Atlantic_Basin", "Pacific_Basin", "Indian_Basin",
                   "Arctic_Basin", "Southern_Ocean_Basin", "Mediterranean_Basin",
                   "Global Ocean", "Global Ocean 65N to 65S",
                   "Global Ocean 15S to 15N"]

    # whether to compute an anomaly with respect to the start of the time series
    computeAnomaly = False

    # Number of points over which to compute moving average(e.g., for monthly
    # output, movingAverageMonths=12 corresponds to a 12-month moving average
    # window)
    movingAverageMonths = 12


    [hovmollerOceanRegionsPotentialTemperature]
    ## options related to plotting time series of temperature vs. depth in ocean
    ## regions

    # colormap
    colormapNameResult = RdYlBu_r
    # whether the colormap is indexed or continuous
    colormapTypeResult = continuous
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': -2., 'vmax': 30.}

    # colormap for differences
    colormapNameDifference = balance
    # whether the colormap is indexed or continuous
    colormapTypeDifference = continuous
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -5., 'vmax': 5.}

    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    contourLevels = 'none'

    # An optional first year for the tick marks on the x axis. Leave commented out
    # to start at the beginning of the time series.

    # firstYearXTicks = 1

    # An optional number of years between tick marks on the x axis.  Leave
    # commented out to determine the distance between ticks automatically.

    # yearStrideXTicks = 1

    # limits on depth, the full range by default
    # yLim = [-6000., 0.]


    [hovmollerOceanRegionsSalinity]
    ## options related to plotting time series of salinity vs. depth in ocean
    ## regions

    # colormap
    colormapNameResult = haline
    # whether the colormap is indexed or continuous
    colormapTypeResult = continuous
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': 30, 'vmax': 39.0}

    # colormap for differences
    colormapNameDifference = balance
    # whether the colormap is indexed or continuous
    colormapTypeDifference = continuous
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.5, 'vmax': 0.5}

    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    contourLevels = 'none'

    # An optional first year for the tick marks on the x axis. Leave commented out
    # to start at the beginning of the time series.

    # firstYearXTicks = 1

    # An optional number of years between tick marks on the x axis.  Leave
    # commented out to determine the distance between ticks automatically.

    # yearStrideXTicks = 1

    # limits on depth, the full range by default
    # yLim = [-6000., 0.]


    [hovmollerOceanRegionsPotentialDensity]
    ## options related to plotting time series of potential density vs. depth in
    ## ocean regions

    # colormap
    colormapNameResult = Spectral_r
    # whether the colormap is indexed or continuous
    colormapTypeResult = continuous
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': 1026.5, 'vmax': 1028.}


    # colormap for differences
    colormapNameDifference = balance
    # whether the colormap is indexed or continuous
    colormapTypeDifference = continuous
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.3, 'vmax': 0.3}

    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    contourLevels = 'none'

    # An optional first year for the tick marks on the x axis. Leave commented out
    # to start at the beginning of the time series.

    # firstYearXTicks = 1

    # An optional number of years between tick marks on the x axis.  Leave
    # commented out to determine the distance between ticks automatically.

    # yearStrideXTicks = 1

    # limits on depth, the full range by default
    # yLim = [-6000., 0.]

The ``[hovmollerOceanRegions]`` section contains a list of ``regionGroups``,
one or more of the :ref:`config_region_groups` defined in
:py:func:`mpas_analysis.shared.regions.compute_region_masks_subtask.get_region_info()`.

For each region group, there is a corresponding section
``[hovmoller<RegionGroup>]``, where ``<RegionGroup>`` is the name of the region
group with spaces removed.  In this section, the  ``fields`` dictionary is used
to specify a list of 3D MPAS fields to average and plot.  The key ``prefix`` is
a convenient name appended to tasks and file names to describe the field.
``mpas`` is the name of the field in MPAS ``timeSeriesStatsMonthly`` output
files.  The ``units`` are the SI units of the field to include on the plot's
color bar  and ``titleName`` is the name of the field to use in its gallery name
and on the title of the plot.

``regionNames`` is a list of regions from the full the region group or
``regionNames = ['all']`` to indicate that all regions should be used.  For
the available regions, see
`Aggregate Existing Features <http://mpas-dev.github.io/geometric_features/stable/aggregation.html>`_.

If ``computeAnomaly = True``, the plots will be the anomaly with respect to the
beginning of the time series (averaged over ``movingAverageMonths`` months).

Each field has a ``[hovmollerOceanRegions<FieldName>]`` section, where
``<FieldName>`` is the associated ``prefix`` but starting with a capital letter.
Each of these sections has a ``yLim`` option that can specify the desired depth
range.  The default is the full range.

For more details on the remaining config options, see
 * :ref:`config_region_groups`
 * :ref:`config_colormaps`
 * :ref:`config_moving_average`
 * :ref:`config_time_axis_ticks`

Example Result
--------------

.. image:: examples/hovmoller_weddell.png
   :width: 500 px
   :align: center
