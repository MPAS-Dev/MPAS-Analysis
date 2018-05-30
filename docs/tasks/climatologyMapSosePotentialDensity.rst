.. _task_climatologyMapSosePotentialDensity:

climatologyMapSosePotentialDensity
==================================

An analysis task for comparing potential density at various depths against
results from the `Southern Ocean State Estimate (SOSE)`_.

Component and Tags::

  component: ocean
  tags: climatology, horizontalMap, sose, potentialDensity, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

    [climatologyMapSosePotentialDensity]
    ## options related to plotting climatology maps of Antarctic
    ## potential density at various levels, including the sea floor against
    ## reference model results and SOSE reanalysis data

    # comparison grid(s)
    # only the Antarctic really makes sense but lat-lon could technically work.
    comparisonGrids = ['antarctic']

    # Times for comparison times (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep,
    # Oct, Nov, Dec, JFM, AMJ, JAS, OND, ANN)
    seasons =  ['ANN','JFM','JAS']

    # list of depths in meters (positive up) at which to analyze, 'top' for the
    # sea surface, 'bot' for the sea floor
    depths = ['top', -200, -400, -600, -800, 'bot']

    # colormap for model/observations
    colormapNameResult = Spectral_r
    # the type of norm used in the colormap
    normTypeResult = linear
    # A dictionary with keywords for the norm
    normArgsResult = {'vmin': 1026.5, 'vmax': 1028.}
    # place the ticks automatically by default
    # colorbarTicksResult = numpy.linspace(1026., 1028., 9)

    # colormap for differences
    colormapNameDifference = balance
    # the type of norm used in the colormap
    normTypeDifference = linear
    # A dictionary with keywords for the norm
    normArgsDifference = {'vmin': -0.3, 'vmax': 0.3}
    # place the ticks automatically by default
    # colorbarTicksDifference = numpy.linspace(-0.3, 0.3, 9)

For more details, see:
 * :ref:`config_colormaps`
 * :ref:`config_seasons`
 * :ref:`config_comparison_grids`

The option ``depths`` is a list of (approximate) depths at which to sample
the potential density field.  A value of ``'top'`` indicates the sea
surface (or the ice-ocean interface under ice shelves) while a value of
``'bot'`` indicates the seafloor.

State Estimate
--------------

:ref:`sose`

Example Result
--------------

.. image:: examples/clim_sose_temp.png
   :width: 720 px
   :align: center

.. _`Southern Ocean State Estimate (SOSE)`: http://sose.ucsd.edu/sose_stateestimation_data_05to10.html
