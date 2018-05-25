.. _task_climatologyMapAntarcticMelt:

climatologyMapAntarcticMelt
===========================

An analysis task for comparison of Antarctic maps of melt rates against
observations from `Rignot et al. (2013)`_.

Component and Tags::

  component: ocean
  tags: climatology, horizontalMap, meltRate, landIceCavities


Configuration Options
---------------------

The following configuration options are available for this task::

  [climatologyMapAntarcticMelt]
  ## options related to plotting horizontally regridded maps of Antarctic
  ## sub-ice-shelf melt rates against reference model results and observations

  # comparison grid(s)
  # only the Antarctic really makes sense but lat-lon could technically work.
  comparisonGrids = ['antarctic']

  # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep,
  # Oct, Nov, Dec, JFM, AMJ, JAS, OND, ANN)
  seasons =  ['JFM', 'JAS', 'ANN']

  # colormap for model/observations
  colormapNameResult = erdc_iceFire_H
  # the type of norm used in the colormap
  normTypeResult = symLog
  # A dictionary with keywords for the norm
  normArgsResult = {'linthresh': 1., 'linscale': 0.5, 'vmin': -100.,
                    'vmax': 100.}
  colorbarTicksResult = [-100., -50., -20., -10., -5., -2., -1., 0., 1., 2.,
                         5., 10., 20., 50., 100.]

  # colormap for differences
  colormapNameDifference = balance
  # the type of norm used in the colormap
  normTypeDifference = symLog
  # A dictionary with keywords for the norm
  normArgsDifference = {'linthresh': 1., 'linscale': 0.5, 'vmin': -100.,
                        'vmax': 100.}
  colorbarTicksDifference = [-100., -50., -20., -10., -5., -2., -1., 0., 1.,
                             2., 5., 10., 20., 50., 100.]

For more details, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`

Observations
------------

:ref:`rignot_melt`

Example Result
--------------

.. image:: examples/ant_melt.png
   :width: 720 px
   :align: center

.. _`Rignot et al. (2013)`: http://doi.org/10.1126/science.1235798