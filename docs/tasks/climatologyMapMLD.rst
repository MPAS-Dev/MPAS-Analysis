.. _task_climatologyMapMLD:

climatologyMapMLD
=================

An analysis task for comparison of global maps of mixed layer depth (MLD)
against observations.

Component and Tags::

  component: ocean
  tags: climatology, horizontalMap, mld, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapMLD]
  ## options related to plotting horizontally remapped climatologies of
  ## mixed layer depth (MLD) against control model results and observations

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

  # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep, Oct,
  # Nov, Dec, JFM, AMJ, JAS, OND, ANN)
  seasons =  ['JFM', 'JAS', 'ANN']

  # comparison grid(s) ('latlon', 'antarctic') on which to plot analysis
  comparisonGrids = ['latlon']

For more details, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`

Observations
------------

:ref:`argo_mld`

Example Result
--------------

.. image:: examples/mld.png
   :width: 500 px
   :align: center
