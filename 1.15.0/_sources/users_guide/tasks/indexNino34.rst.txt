.. |n~| unicode:: U+00F1
   :trim:

.. _task_indexNino34:

indexNino34
===========

An analysis task for plotting both time series and spectra of the El Ni |n~| o
3.4 Climate Index against observations.

Component and Tags::

  component: ocean
  tags: timeSeries, index, nino, publicObs

Configuration Options
---------------------

The following configuration options are available for this task::

  [indexNino34]
  ## options related to plotting time series of the El Nino 3.4 index

  # Specified region for the Nino Index,'nino3', 'nino4', or 'nino3.4'
  # The indexNino34 routine only accepts one value at a time
  region = nino3.4

  # Data source to read for comparison.  There are two options
  # 1 - ERS_SSTv4 -- Extended Reconstructed Sea Surface Temperature -- 1854 - 2016
  # 2 - HADIsst -- Hadley Center analysis -- 1870 - 2016
  observationData = HADIsst

While the default is the El Ni |n~| o 3.4 region, you may select among
``nino3``, ``nino4``, and ``nino3.4``.  See :ref:`config_regions` for more
information about regions in MPAS-Analyis.

By default, observations are taken from the Hadley Center analysis.  To use
the Extended Reconstructed Sea Surface Temperature (ERS SSTv4), set
``observationData = ERS_SSTv4``.

Observations
------------

* :ref:`hadisst_nino`
* :ref:`ers_sst_nino`


Example Result
--------------

.. image:: examples/nino.png
   :width: 500 px
   :align: center
