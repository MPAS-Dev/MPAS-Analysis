.. _task_climatologyMapSeaIceSnowDepthSH:

climatologyMapSeaIceSnowDepthSH
=======================================

An analysis task for plotting maps of Antarctic sea ice cell-averaged snow depth.

Component and Tags::

  component: seaIce
  tags: climatology, horizontalMap, seaIceSnowDepth, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapSeaIceSnowDepthSH]
  ## options related to plotting horizontally remapped climatologies of
  ## snow depth against control model results and observations
  ## in the southern hemisphere (SH)

  # colormap for model/observations
  colormapNameResult = devon
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [20, 80, 110, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 1, 1.5]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 32, 64, 96, 112, 128, 128, 144, 160, 192, 224, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-1.5, -1, -0.5, -0.25, -0.05, 0, 0.05, 0.25, 0.5, 1, 1.5]

  # Times for comparison
  seasons =  ['FM', 'ON']

  # comparison grid(s) (typically 'arctic_extended') on which to plot analysis
  comparisonGrids = ['antarctic_extended']

  # arrange subplots vertically?
  vertical = False

The option ``vertical = True`` can be used to plot 3 panels one above another
(resulting in a tall, thin image) rather than next to each other, the default
(resulting in a short, wide image).

For details on the remaining configration options, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`

Example Result
--------------

.. image:: examples/ice_snowdepth_sh.png
   :width: 720 px
   :align: center
