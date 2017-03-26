#############
API reference
#############

This page provides an auto-generated summary of the MPAS-Analysis API. For
more details and examples, refer to the relevant chapters in the main part of
the documentation.

Top-level script: run_analysis.py
=================================
.. autosummary::
   :toctree: generated/

   run_analysis.check_path_exists
   run_analysis.check_generate
   run_analysis.analysis


Analysis tasks
==============

Ocean tasks
-----------

.. currentmodule:: mpas_analysis.ocean

.. autosummary::
   :toctree: generated/

   meridional_overturning_circulation.moc_streamfunction
   nino34_index.nino34_index
   ocean_modelvsobs.ocn_modelvsobs
   ohc_timeseries.ohc_timeseries
   sst_timeseries.sst_timeseries

Sea ice tasks
-------------

.. currentmodule:: mpas_analysis.sea_ice

.. autosummary::
   :toctree: generated/

   modelvsobs.seaice_modelvsobs
   timeseries.seaice_timeseries


Configuration
==============
.. currentmodule:: mpas_analysis.configuration.MpasAnalysisConfigParser

.. autosummary::
   :toctree: generated/

   MpasAnalysisConfigParser.getWithDefault
   MpasAnalysisConfigParser.getExpression


Shared modules
==============

Reading MPAS Datasets
---------------------
.. currentmodule:: mpas_analysis.shared.mpas_xarray

.. autosummary::
   :toctree: generated/

   mpas_xarray.open_multifile_dataset
   mpas_xarray.preprocess
   mpas_xarray.remove_repeated_time_index
   mpas_xarray.subset_variables

.. currentmodule:: mpas_analysis.shared.generalized_reader

.. autosummary::
   :toctree: generated/

   generalized_reader.open_multifile_dataset


Climatology
-----------
.. currentmodule:: mpas_analysis.shared.climatology

.. autosummary::
   :toctree: generated/

   climatology.compute_monthly_climatology
   climatology.compute_seasonal_climatology
   climatology.get_mpas_climatology_file_names
   climatology.get_observation_climatology_file_names
   climatology.write_mpas_mapping_file
   climatology.write_observations_mapping_file

Interpolation
-------------
.. currentmodule:: mpas_analysis.shared.interpolation

.. autosummary::
   :toctree: generated/

   interpolate.build_remap_weights
   interpolate.remap
   scrip.mpas_file_to_scrip
   scrip.lat_lon_file_to_scrip
   scrip.lat_lon_array_to_scrip

Namelist and Streams Files
--------------------------
.. currentmodule:: mpas_analysis.shared.io.namelist_streams_interface

.. autosummary::
   :toctree: generated/

   convert_namelist_to_dict
   NameList.__init__
   NameList.__getattr__
   NameList.__getitem__
   NameList.get
   NameList.getint
   NameList.getfloat
   NameList.getbool

   StreamsFile.__init__
   StreamsFile.read
   StreamsFile.readpath
   StreamsFile.has_stream
   StreamsFile.find_stream

I/O Utilities
-------------
.. currentmodule:: mpas_analysis.shared

.. autosummary::
   :toctree: generated/

   io.utility.paths
   io.utility.make_directories
   io.utility.build_config_full_path

Plotting
--------
.. currentmodule:: mpas_analysis.shared.plot

.. autosummary::
   :toctree: generated/

   plotting.nino34_spectra_plot
   plotting.nino34_timeseries_plot
   plotting.timeseries_analysis_plot
   plotting.timeseries_analysis_plot_polar
   plotting.plot_polar_comparison
   plotting.plot_global_comparison
   plotting.plot_vertical_section
   plotting.setup_colormap

Timekeeping
-----------
.. currentmodule:: mpas_analysis.shared.timekeeping

.. autosummary::
   :toctree: generated/

   utility.get_simulation_start_time
   utility.string_to_datetime
   utility.string_to_relative_delta
   utility.string_to_days_since_date
   utility.days_to_datetime
   utility.datetime_to_days
   utility.date_to_days
   MpasRelativeDelta.MpasRelativeDelta

