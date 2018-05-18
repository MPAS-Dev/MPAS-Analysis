.. _config_seasons:

Seasons
=======

Nearly all analysis tasks that produce climatology plots include a
configuration option for specifying a list of seasons to plot::

  # Months or seasons to plot (Jan, Feb, Mar, Apr, May, Jun, Jul, Aug, Sep,
  # Oct, Nov, Dec, JFM, AMJ, JAS, OND, ANN)
  seasons =  ['JFM', 'JAS', 'ANN']

Valid seasons include the three-letter abbreviations of each month (e.g.
``Jan``), several 2- and 3-month seasons specified by consecutive first letters
of each month name (``JFM``, ``AMJ``, ``JAS``, ``OND``, ``ON``, ``FM``, ``DJF``
and ``JJA``), and ``ANN`` for all 12 months.

If seasons other than these are needed, please post an issue on `GitHub`_ or
contact the developers.

.. _`GitHub`: https://github.com/MPAS-Dev/MPAS-Analysis/issues