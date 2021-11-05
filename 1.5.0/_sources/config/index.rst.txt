.. |n~| unicode:: U+00F1
   :trim:

.. _config_index:

Index
=====

The ``[index]`` section of a configuration file contains options used to
determine the start and end years of climate indices (such as El Ni |n~| o
3.4)::

  [index]
  ## options related to producing nino index.

  # start and end years for El Nino 3.4 analysis. Use endYear = end to indicate
  # that the full range of the data should be used.  If errorOnMissing = False,
  # the start and end year will be clipped to the valid range.  Otherwise, out
  # of bounds values will lead to an error.  In a "control" config file used in
  # a "main vs. control" analysis run, the range of years must be valid and
  # cannot include "end" because the original data may not be available.
  startYear = 1
  endYear = end

Start and End Year
------------------

A custom config file should specify a start and end year for time axis.
If ``errorOnMissing = False`` in the ``input`` section and the start or end
year is beyond the range of the simulation, the range will be reduced to those
dates with available data and a warning message will be displayed.  If
``errorOnMissing = True``, out of range year will produce an error.


