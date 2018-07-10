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

  # start and end years for the nino 3.4 analysis.  Using out-of-bounds values
  #   like start_year = 1 and end_year = 9999 will be clipped to the valid range
  #   of years, and is a good way of insuring that all values are used.
  # For valid statistics, index times should include at least 30 years
  startYear = 1
  endYear = 9999

Start and End Year
------------------

A custom config file should specify a start and end year for climate indices.
If the end year is beyond the range of the simulation (typically true for the
default ``endYear = 9999``), the range will be reduced to those dates with
available data and a warning message will be displayed.


