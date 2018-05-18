.. _config_plot:

Plot
====

The ``[plot]`` section of a configuration file contains options that define
default properties for all plots produce by MPAS-Analysis::

  [plot]
  ## options related to plotting that are the defaults across all analysis
  ## modules

  # set to true if you want plots to be displayed (one by one) to the screen in
  # addition to being written out to png files
  # Note: displayToScreen = True seems to hang on Edison on large data sets,
  # so suggested use is just for debugging either locally or with small data sets
  displayToScreen = False

  # font size on axes
  axisFontSize = 16
  # title font properties
  titleFontSize = 20
  titleFontColor = black
  titleFontWeight = normal

  # the dots per inch of output figures
  dpi = 200


The option ``displayToScreen = True`` was used in the early stages of
development for debugging and has not been tested recently (meaning it likely
no longer works).  In theory, this would allow all plot to be displayed to
the screen.  This functionality almost certainly does not work when running
with parallel tasks (see :ref:`config_execute`).

The options for title fint size, color and weight as well as axis font size
specify properties of these parts of each plot.  Sizes are given in points.

The value of ``dpi`` specifies the resolution of the images written out by
MPAS-Analysis (in dots per inch).  The default produces large images that
are appropriate for zooming in substantially and may be sufficient for
publication.  They are large (but not entirely unmanageable) for the web.
