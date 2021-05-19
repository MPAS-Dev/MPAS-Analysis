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

  # the default font size for axis ticks, etc.
  defaultFontSize = 10

  # title and axis font properties for single-panel plots
  titleFontSize = 16
  titleFontColor = black
  titleFontWeight = normal
  axisFontSize = 12

  # font size for cartopy grid labels
  cartopyGridFontSize = 12

  # the dots per inch of output figures
  dpi = 200

  # Write out PDFs in addition to PNGs?
  pdf = False

The options for title font size, color and weight as well as axis font size
specify properties of these parts of each plot.  The default font size covers
the axis tick marks, color-bar ticks and axis labels, contour labels, etc.
The cartopy grid font size are for the labels (either interior or along the
plot boundary) for cartopy labels. Sizes are given in points.

The value of ``dpi`` specifies the resolution of the images written out by
MPAS-Analysis (in dots per inch).  The default produces large images that
are appropriate for zooming in substantially and may be sufficient for
publication.  They are large (but not entirely unmanageable) for the web.

You can set ``pdf = True`` to write out PDF files in the plots subdirectory
along with PNG files.  The PDFs are not copied to the HTML folder.

Many types of individual plots, including climatologies, transects, Hovmoller
plots, and most time series, also support setting the ``defaultFontSize`` and
``titleFontSize`` config options just of that type of plot.
