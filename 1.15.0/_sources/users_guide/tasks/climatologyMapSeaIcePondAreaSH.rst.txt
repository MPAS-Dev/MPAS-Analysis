.. _task_climatologyMapSeaIcePondAreaSH:

climatologyMapSeaIcePondAreaSH
=======================================

An analysis task for plotting maps of Antarctic sea ice pond area fraction.

Component and Tags::

  component: seaIce
  tags: climatology, horizontalMap, seaIcePondArea, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapSeaIcePondAreaSH]
  # options related to plotting horizontally remapped climatologies of
  # snow melt against control model results and observations
  # in the southern hemisphere (SH)

  # colormap for model/observations
  colormapNameResult = deep
  # whether the colormap is indexed or continuous
  colormapTypeResult = indexed
  # color indices into colormapName for filled contours
  colormapIndicesResult = [20, 80, 110, 140, 170, 200, 230, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsResult = [0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4]

  # colormap for differences
  colormapNameDifference = balance
  # whether the colormap is indexed or continuous
  colormapTypeDifference = indexed
  # color indices into colormapName for filled contours
  colormapIndicesDifference = [0, 32, 64, 96, 128, 128, 160, 192, 224, 255]
  # colormap levels/values for contour boundaries
  colorbarLevelsDifference = [-0.2, -0.1, -0.05, -0.01, -0.005, 0, 0.005, 0.01, 0.05, 0.1, 0.2]

  # Months or seasons to plot
  seasons =  ['FM', 'ON']

  # comparison grid(s) (typically 'antarctic_extended') on which to plot analysis
  comparisonGrids = ['antarctic_extended']

  # a list of prefixes describing the sources of the observations to be used
  observationPrefixes = ['']

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

.. image:: examples/ice_pondarea_sh.png
   :width: 720 px
   :align: center
