.. _task_climatologyMapWoa:

climatologyMapWoa
=================

An analysis task for comparing potential temperature and salinity
at various depths against WOA18 climatology.

Component and Tags::

  component: ocean
  tags: climatology, horizontalMap, woa, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapWoa]
  ## options related to plotting climatology maps of Temperature and Salinity
  ## fields at various levels, including the sea floor against control model
  ## results and WOA climatological data

  # comparison grid(s) ('latlon', 'antarctic', 'arctic') on which to plot analysis
  comparisonGrids = ['arctic']

  # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
  # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
  seasons =  ['ANN','JFM','JAS']

  # list of depths in meters (positive up) at which to analyze, 'top' for the
  # sea surface. Note that, for seasons='ANN', depths can be as deep as 5500 m,
  # otherwise use a maximum depth of 1500 m.
  depths = ['top', -500, -1000]

  # a list of fields top plot for each transect.  All supported fields are listed
  # below
  fieldList = ['temperature', 'salinity']

  [climatologyMapWoaTemperature]
  ## options related to plotting climatology maps of potential temperature
  ## at various levels, including the sea floor against control model results
  ## and WOA18 climatological data

  # colormap for model/observations
  colormapNameResult = RdYlBu_r
  # whether the colormap is indexed or continuous
  colormapTypeResult = continuous
  # the type of norm used in the colormap
  normTypeResult = linear
  # A dictionary with keywords for the norm
  normArgsResult = {'vmin': -2., 'vmax': 2.}
  # place the ticks automatically by default
  # colorbarTicksResult = numpy.linspace(-2., 2., 9)

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = continuous
  # the type of norm used in the colormap
  normTypeDifference = linear
  # A dictionary with keywords for the norm
  normArgsDifference = {'vmin': -2., 'vmax': 2.}
  # place the ticks automatically by default
  # colorbarTicksDifference = numpy.linspace(-2., 2., 9)

  [climatologyMapWoaSalinity]
  ## options related to plotting climatology maps of salinity
  ## at various levels, including the sea floor against control model results 
  ## and WOA18 climatological data

  # colormap for model/observations
  colormapNameResult = haline
  # whether the colormap is indexed or continuous
  colormapTypeResult = continuous
  # the type of norm used in the colormap
  normTypeResult = linear
  # A dictionary with keywords for the norm
  normArgsResult = {'vmin': 33.8, 'vmax': 35.0}
  # place the ticks automatically by default
  # colorbarTicksResult = numpy.linspace(34.2, 35.2, 9)

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = continuous
  # the type of norm used in the colormap
  normTypeDifference = linear
  # A dictionary with keywords for the norm
  normArgsDifference = {'vmin': -0.5, 'vmax': 0.5}
  # place the ticks automatically by default
  # colorbarTicksDifference = numpy.linspace(-0.5, 0.5, 9)

For more details, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`

The option ``depths`` is a list of (approximate) depths at which to sample
the temperature and salinity fields.  A value of ``'top'`` indicates the sea
surface. Note that, for the annual climatology, WOA18 data is available down
to 5500 m, whereas, for the seasonal or monthly climatologies, WOA18 data
is only available down to 1500 m.

Observations
------------

:ref:`woa18_t_s`

Example Result
--------------

.. image:: examples/clim_woa_temp.png
   :width: 500 px
   :align: center
