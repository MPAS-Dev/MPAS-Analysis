#!/usr/bin/env python

"""
Runs MPAS-Analysis via configuration file `config.analysis` specifying analysis
options.

Author: Xylar Asay-Davis, Phillip J. Wolfram
Last Modified: 9/15/2016
"""

import os
import sys
import matplotlib as mpl

from mpas_analysis.configuration.MpasAnalysisConfigParser import MpasAnalysisConfigParser

def path_existence(config, section, option, ignorestr=None): #{{{
    inpath = config.get(section, option)
    if not (os.path.isdir(inpath) or os.path.isfile(inpath)):
        # assumes that path locations of ignorestr won't return an error, e.g.,
        # ignorestr="none" is a key word to indicate the path or file is
        # optional and is not needed
        if inpath == ignorestr:
            return False
        errmsg = "Path %s not found. Exiting..." % inpath
        raise SystemExit(errmsg)
    return inpath #}}}

def analysis(config): #{{{

    # Checks on directory/files existence:
    indir = path_existence(config, 'paths', 'archive_dir_ocn')

    if config.get('case', 'ref_casename_v0') != 'None':
        path_existence(config, 'paths','ref_archive_v0_ocndir')
        path_existence(config, 'paths','ref_archive_v0_seaicedir')

    use_seaice = config.getboolean('seaice_timeseries', 'generate')
    seaice_compare_obs = config.getboolean('seaice_timeseries', 'compare_with_obs')
    seaice_modelvsobs = config.getboolean('seaice_modelvsobs', 'generate')
    if (use_seaice and seaice_compare_obs) or seaice_modelvsobs:
        # we will need sea-ice observations.  Make sure they're there
        for obsfile in ['obs_iceareaNH', 'obs_iceareaSH', 'obs_icevolNH', 'obs_icevolSH']:
            path_existence('seaIceData', obsfile, ignorestr='none')

    # choose the right rendering backend, depending on whether we're displaying
    # to the screen
    if not config.getboolean('plot', 'displayToScreen'):
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    # analysis can only be imported after the right MPL renderer is selected

    #GENERATE OCEAN DIAGNOSTICS
    if config.getboolean('ohc_timeseries','generate'):
        print ""
        print "Plotting OHC time series..."
        from mpas_analysis.ocean.ohc_timeseries import ohc_timeseries
        ohc_timeseries(config)

    if config.getboolean('sst_timeseries','generate'):
        print ""
        print "Plotting SST time series..."
        from mpas_analysis.ocean.sst_timeseries import sst_timeseries
        sst_timeseries(config)

    if config.getboolean('nino34_timeseries','generate'):
        print ""
        print "Plotting Nino3.4 time series..."
        #from mpas_analysis.ocean.nino34_timeseries import nino34_timeseries
        #nino34_timeseries(config)

    if config.getboolean('mht_timeseries','generate'):
        print ""
        print "Plotting Meridional Heat Transport (MHT)..."
        #from mpas_analysis.ocean.mht_timeseries import mht_timeseries
        #mht_timeseries(config)


    if config.getboolean('moc_timeseries','generate'):
        print ""
        print "Plotting Meridional Overturning Circulation (MOC)..."
        #from mpas_analysis.ocean.moc_timeseries import moc_timeseries
        #moc_timeseries(config)

    if config.getboolean('sst_modelvsobs','generate'):
        print ""
        print "Plotting 2-d maps of SST climatologies..."
        from mpas_analysis.ocean.ocean_modelvsobs import sst_modelvsobs
        sst_modelvsobs(config)
        
    if config.getboolean('mld_modelvsobs','generate'):
        print ""
        print "Plotting 2-d maps of MLD climatologies..."
        from mpas_analysis.ocean.ocean_modelvsobs import mld_modelvsobs
        mld_modelvsobs(config)


    #GENERATE SEA-ICE DIAGNOSTICS
    if config.getboolean('seaice_timeseries','generate'):
        print ""
        print "Plotting sea-ice area and volume time series..."
        from mpas_analysis.sea_ice.timeseries import seaice_timeseries
        seaice_timeseries(config)

    if config.getboolean('seaice_modelvsobs','generate'):
        print ""
        print "Plotting 2-d maps of sea-ice concentration and thickness climatologies..."
        from mpas_analysis.sea_ice.modelvsobs import seaice_modelvsobs
        seaice_modelvsobs(config)


    #GENERATE LAND-ICE DIAGNOSTICS


    if config.getboolean('plot','displayToScreen'):
        plt.show()


    return #}}}

if __name__ == "__main__":

    # process command line arguments and run analysis from configuration
    if len(sys.argv) <= 1:
        print "usage: %s <in_config_file> [<in_config_file2>]" % sys.argv[0]
        exit(1)

    configFileNames = sys.argv[1:]
    config = MpasAnalysisConfigParser()
    config.read(configFileNames)

    analysis(config)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
