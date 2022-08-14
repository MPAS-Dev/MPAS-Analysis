.. _config_time_axis_ticks:

Time-Axis Tick Marks
====================

By default tick marks on the time (x) axis are determined automatically by
``matplotlib``.  You can explicitly control them by setting a first year and
a stride (skip) in years.  For example::

  firstYearXTicks = 2
  yearStrideXTicks = 2

will place the first tick mark at simulation year 2 and will give a tick mark
every 2 years.


