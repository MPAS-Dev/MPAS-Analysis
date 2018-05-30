.. _config_plot:

Plot
====

The ``[plot]`` section of a configuration file contains options that define
default properties for all plots produce by MPAS-Analysis::

  [plot]
  ## options related to plotting that are the defaults across all analysis
  ## modules

  # font size on axes
  axisFontSize = 16
  # title font properties
  titleFontSize = 20
  titleFontColor = black
  titleFontWeight = normal

  # the dots per inch of output figures
  dpi = 200

The options for title fint size, color and weight as well as axis font size
specify properties of these parts of each plot.  Sizes are given in points.

The value of ``dpi`` specifies the resolution of the images written out by
MPAS-Analysis (in dots per inch).  The default produces large images that
are appropriate for zooming in substantially and may be sufficient for
publication.  They are large (but not entirely unmanageable) for the web.
