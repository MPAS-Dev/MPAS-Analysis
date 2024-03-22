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

  # start and end years for El Nino 3.4 analysis.  Out-of-bounds values will lead
  # to an error.
  startYear = 1
  endYear = 20

Start and End Year
------------------

A custom config file should specify a start and end year for time axis.
Out of range year will produce an error.


