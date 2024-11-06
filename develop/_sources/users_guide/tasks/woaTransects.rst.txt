.. |deg| unicode:: U+00B0 .. degree sign
   :trim:

.. _task_woaTransects:

woaTransects
============

An analysis task for computing meridional transects of MPAS fields at evenly
spaced latitudes around Antarctica and comparing them with results from the
`World Ocean Atlas 2023 (WOA23)`_.

Component and Tags::

    component: ocean
    tags: climatology, transect, woa, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

    [woaTransects]
    ## options related to plotting model vs. World Ocean Atlas 2023 (WOA23)
    ## transects.

    # Times for comparison times (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
    # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
    seasons =  ['ANN', 'JFM', 'JAS']

    # The approximate horizontal resolution (in km) of each transect.  Latitude/
    # longitude between observation points will be subsampled at this interval.
    # Use 'obs' to indicate no subsampling. Use 'mpas' to indicate plotting of
    # model data on the native grid, in which case comparison with observations
    # will take place on the observation grid.
    #horizontalResolution = mpas
    #horizontalResolution = obs
    horizontalResolution = 5

    # The name of the vertical comparison grid.  Valid values are 'mpas' for the
    # MPAS vertical grid, 'obs' to use the locations of observations or
    # any other name if the vertical grid is defined by 'verticalComparisonGrid'
    #verticalComparisonGridName = mpas
    #verticalComparisonGridName = obs
    verticalComparisonGridName = uniform_10_to_1500m_at_10m

    # The vertical comparison grid if 'verticalComparisonGridName' is not 'mpas' or
    # 'obs'.  This should be numpy array of (typically negative) elevations (in m).
    # The first and last entries are used as axis bounds for 'mpas' and 'obs'
    # vertical comparison grids
    verticalComparisonGrid = numpy.linspace(-10, -1500, 150)

    # A range for the y axis (if any)
    verticalBounds = []

    # The minimum weight of a destination cell after remapping. Any cell with
    # weights lower than this threshold will therefore be masked out.
    renormalizationThreshold = 0.01

    # min and max latitude of transects
    minLat = -80
    maxLat = -60

    # longitudes of transects
    # Default transects are at Filchner, Bellingshausen Sea, Thwaites, Ross, Totten,
    # Amery, Fimbul
    longitudes = [318., 280., 253., 187., 117., 75., 0.]

    # a list of fields top plot for each transect.  All supported fields are listed
    # below.
    fieldList = ['temperature', 'salinity']


    [woaTemperatureTransects]
    ## options related to plotting WOA23 transects of potential temperature

    # colormap for model/observations
    colormapNameResult = RdYlBu_r
    # whether the colormap is indexed or continuous
    colormapTypeResult = continuous
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': -2.0, 'vmax': 2.0}
    # color indices into colormapName for filled contours
    #colormapIndicesResult = [0, 40, 80, 110, 140, 170, 200, 230, 255]
    # colormap levels/values for contour boundaries
    #colorbarLevelsResult = [0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5, 6]
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(0.0, 6.0, 9)
    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    #contourLevelsResult = np.arange(0.5, 6.0, 1.0)
    contourLevelsResult = 'none'

    # colormap for differences
    colormapNameDifference = balance
    # whether the colormap is indexed or continuous
    colormapTypeDifference = continuous
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -2.0, 'vmax': 2.0}
    # color indices into colormapName for filled contours
    #colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
    # colormap levels/values for contour boundaries
    #colorbarLevelsDifference = [-2, -1.5, -1.25, -1, -0.2, 0, 0.2, 1, 1.25, 1.5, 2]
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-2.0, 2.0, 9)
    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    #contourLevelsDifference = np.arange(-1.8, 2.0, 0.4)
    contourLevelsDifference = 'none'


    [woaSalinityTransects]
    ## options related to plotting WOA23 transects of salinity

    # colormap for model/observations
    colormapNameResult = haline
    # whether the colormap is indexed or continuous
    colormapTypeResult = continuous
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': 33.8, 'vmax': 35.0}
    # color indices into colormapName for filled contours
    #colormapIndicesResult = [0, 40, 80, 110, 140, 170, 200, 230, 255]
    # colormap levels/values for contour boundaries
    #colorbarLevelsResult = [34, 34.3, 34.5, 34.65, 34.675, 34.7, 34.725, 34.75, 34.8, 35]
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(34.0, 35.0, 9)
    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    #contourLevelsResult = np.arange(34.1, 35.0, 0.1)
    contourLevelsResult = 'none'

    # colormap for differences
    colormapNameDifference = balance
    # whether the colormap is indexed or continuous
    colormapTypeDifference = continuous
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.5, 'vmax': 0.5}
    # color indices into colormapName for filled contours
    #colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
    # colormap levels/values for contour boundaries
    #colorbarLevelsDifference = [-0.5, -0.2, -0.1, -0.05, -0.02, 0,  0.02, 0.05, 0.1, 0.2, 0.5]
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-0.5, 0.5, 9)
    # contour line levels (use [] for automatic contour selection, 'none' for no
    # contour lines)
    #contourLevelsDifference = numpy.linspace(-0.6, 0.6, 9)
    contourLevelsDifference = 'none'


The options ``minLat`` and ``maxLat`` determine the start and end of each
meridional transect (in degrees).  The option ``longitudes`` is a list or
numpy array of longitudes for each transect, e.g.::

  longitudes = numpy.linspace(0, 330, 12)

produces 12 transects spaced every 30 |deg|.

The user can select only to plot a subset of the supported fields by adding
only the desired field names to ``fieldList``.  Only potential temperature
and salinity are currently available.

Ater the ``woaTransects`` section, there is a section for each supported field
specifying the information related to the colormap.

For details on remaining configuration options, see:
 * :ref:`config_transects`
 * :ref:`config_remapping`
 * :ref:`config_colormaps`
 * :ref:`config_seasons`

Observations
------------

:ref:`woa23_t_s`

Example Result
--------------

.. image:: examples/woa23_transect.png
   :width: 500 px
   :align: center

.. _`World Ocean Atlas 2023 (WOA23)`: https://www.ncei.noaa.gov/products/world-ocean-atlas
