.. _config_diagnostics:

Diagnostics
===========

The ``[diagnostics]`` section of a configuration file contains options related
to paths containing observations, region-mask files and mapping files used to
interpolate MPAS data and observations to common reference grids::

  [diagnostics]
  ## config options related to observations, mapping files and region files used
  ## by MPAS-Analysis in diagnostics computations.

  # The base path to the diagnostics directory.  Typically, this will be a shared
  # directory on each E3SM supported machine (see the example config files for
  # its location).  For other machines, this would be the directory pointed to
  # when running "download_analysis_data.py" to get the public observations,
  # mapping files and region files.
  baseDirectory = /path/to/diagnostics

  # A second directory where custom diagonstics data such as mapping files and
  # regions masks for unsupported grids can be found.  The structure of
  # subdirectories in this directory must be the same as baseDirectory
  customDirectory = none

  # Directory for mapping files (if they have been generated already). If mapping
  # files needed by the analysis are not found here, they will be generated and
  # placed in the output mappingSubdirectory.  The user can supply an absolute
  # path here to point to a path that is not within the baseDirectory above.
  mappingSubdirectory = mpas_analysis/maps

  # Directory for region mask files. The user can supply an absolute path here to
  # point to a path that is not within the baseDirectory above.
  regionMaskSubdirectory = mpas_analysis/region_masks

Diagnostics Directories
-----------------------

The ``baseDirectory`` is the location where files were downloaded with the
``download_analysis_data.py``.  If the user is on an E3SM supported machine,
this data has already been downloaded to a shared location (see example config
files in the subdirectories of the ``configs`` directory in the MPAS-Analysis
repository).

The ``customDirectory``, if it is not ``none`` is another directory where
observations, mapping files and region masks may be stored.  This is useful for
runs on non-standard grids or for testing out new observations.

The remaining options point to the subdirectories for mapping files (see
below) and region masks (see :ref:`config_colormaps`), respectively.
Typically, there is no reason to change ``mappingSubdirectory`` or
``regionMaskSubdirectory``, as these are the standard subdirectories created
when these files are downloaded from the `E3SM public data repository`_.

.. _config_mapping_files:

Mapping Files
-------------

Mapping files are used in many MPAS-Analysis tasks to remap from either the
native MPAS mesh or an observations grid to a comparison grid (see
:ref:`config_comparison_grids`).  By default, these mapping files are generated
on the fly as they are needed.  This can be a time-consuming process,
especially for high resolution meshes, so it is useful to store a cache of
these mapping files for reuse.  Mapping files at three standard resolutions
are avaliable on the `E3SM public data repository`_.  The mapping files for
the two coarser resolution meshes will be downloaded automatically along with
the publicly available observations. (See the :ref:`quick_start` for details
on downloading this data.)

If you notice that MPAS-Analysis is generating mapping files on the fly each
time you run, you may wish to copy them from the mapping files output
directory (the subdirectory ``mapping/`` inside the output base directory) to
your mapping files cache directory.

.. _`E3SM public data repository`: https://web.lcrc.anl.gov/public/e3sm/diagnostics/
