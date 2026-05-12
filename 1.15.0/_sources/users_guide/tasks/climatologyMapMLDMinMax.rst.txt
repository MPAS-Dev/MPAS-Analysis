.. _task_climatologyMapMLDMinMax:

climatologyMapMLDMinMax
=======================

An analysis task for plotting climatologies of monthly min and max of
mixed layer depth (MLD).

Component and Tags::

  component: ocean
  tags: climatology, horizontalMap, mld, publicObs, min, max

Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapMLDMinMax]
  ## options related to plotting horizontally remapped climatologies of
  ## monthly min/max mixed layer depth (MLD), optionally against control model
  ## results

  # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
  # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
  seasons =  ['JFM', 'JAS', 'ANN']

  # comparison grid(s) ('latlon', 'antarctic') on which to plot analysis
  comparisonGrids = ['latlon']


  [climatologyMapDensityMLDMinMax]
  ## color mpas for density-threshold MLD min/max plots in the same figure

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = viridis
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 40, 80, 110, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [0, 1, 5, 10, 20, 30, 50, 80, 110, 150]

  [climatologyMapDensityMLDMin]
  ## color mpas for density MLD min plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]

  [climatologyMapDensityMLDMax]
  ## color mpas for density MLD max plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]


  [climatologyMapTemperatureMLDMinMax]
  ## color mpas for temperature-threshold MLD min/max plots in the same figure

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = viridis
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 40, 80, 110, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [0, 1, 5, 10, 20, 30, 50, 80, 110, 150]

  [climatologyMapTemperatureMLDMin]
  ## color mpas for temperature MLD min plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]

  [climatologyMapTemperatureMLDMax]
  ## color mpas for temperature MLD max plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]


  [climatologyMapBLDMinMax]
  ## color mpas for boundary-layer depth min/max plots in the same figure

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = viridis
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 40, 80, 110, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [0, 1, 5, 10, 20, 30, 50, 80, 110, 150]

  [climatologyMapBLDMin]
  ## color mpas for BLD min plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]

  [climatologyMapBLDMax]
  ## color mpas for BLD max plots in separate figures for main vs. control

  # colormap for model/observations
  colormapNameResult = viridis
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [0, 20, 50, 80, 120, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [10., 30., 50., 75., 100., 150., 200., 400., 800.]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 28, 57, 85, 113, 128, 128, 142, 170, 198, 227, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-150, -80, -30, -10, -1, 0, 1, 10, 30, 80, 150]

Config sections ``climatologyMapDensityMLD*`` apply to plots of the MLD based
on a density threshold, while ``climatologyMapTemperatureMLD*`` are those
related to MLD with a temperature threshold, and ``climatologyMabBLD*`` are of
controling plots of the boundary-layer depth (BLD).  Each of these has a
section ending in ``MinMax`` that controls plots in when not comparing a main
reference to observations, in which case the maximum field is in the top panel,
minimum in the middle panel and the difference in the bottom panel.  When
comparing with a reference run, min. and max. plots are performed separately
along with their corresponding fields from the reference run and the difference
between main and reference.  In these cases, a different color map may be
appropriate (e.g with positive and negative values for the difference, rather
than just positive for max. minus min.).

For more details, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`
