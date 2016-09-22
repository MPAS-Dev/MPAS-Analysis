import os
import sys
import matplotlib as mpl

from mpas_analysis.configuration.MpasAnalysisConfigParser import MpasAnalysisConfigParser

if len(sys.argv) <= 1:
    print "usage: %s <in_config_file> [<in_config_file2>]"%sys.argv[0]
    exit(1)

configFileNames = sys.argv[1:]

config = MpasAnalysisConfigParser()

config.read(configFileNames)

# Checks on directory/files existence:
indir = config.get('paths','archive_dir_ocn')
if not os.path.isdir(indir):
    raise SystemExit("Model directory %s not found. Exiting..." % indir)

ref_casename_v0 = config.get('case','ref_casename_v0')
if ref_casename_v0 != "None":
    # we will need model data.  Make sure it's there
    indir_v0data = config.get('paths','ref_archive_v0_ocndir')
    if not os.path.isdir(indir_v0data):
        raise SystemExit("ref_archive_v0_ocndir directory %s not found. Exiting..." % indir_v0data)
    indir_v0data = config.get('paths','ref_archive_v0_seaicedir')
    if not os.path.isdir(indir_v0data):
        raise SystemExit("ref_archive_v0_seaicedir directory %s not found. Exiting..." % indir_v0data)

if ((config.getboolean('seaice_timeseries', 'generate')
            and config.getboolean('seaice_timeseries', 'compare_with_obs'))
        or config.getboolean('seaice_modelvsobs', 'generate')):
    # we will need sea-ice observations.  Make sure they're there
    for obsFile in ['obs_iceareaNH', 'obs_iceareaSH', 'obs_icevolNH', 'obs_icevolSH']:
        obs_filename = config.get('seaIceData',obsFile)
        if (obs_filename != 'none') and not os.path.isfile(obs_filename):
            raise SystemExit("Obs file %s not found. Exiting..." % obs_filename)

# choose the right rendering backend, depending on whether we're displaying
# to the screen
if not config.getboolean('plot','displayToScreen'):
    mpl.use('Agg')
import matplotlib.pyplot as plt

# these only get imported after we have the right MPL renderer selected
from mpas_analysis.ocean.ohc_timeseries import ohc_timeseries
from mpas_analysis.ocean.sst_timeseries import sst_timeseries
#from mpas_analysis.ocean.nino34_timeseries import nino34_timeseries
#from mpas_analysis.ocean.mht_timeseries import mht_timeseries
#from mpas_analysis.ocean.moc_timeseries import moc_timeseries
from mpas_analysis.ocean.sst_modelvsobs import sst_modelvsobs

from mpas_analysis.sea_ice.timeseries import seaice_timeseries
from mpas_analysis.sea_ice.modelvsobs import seaice_modelvsobs


#GENERATE OCEAN DIAGNOSTICS
if config.getboolean('ohc_timeseries','generate'):
    print ""
    print "Plotting OHC time series..."
    ohc_timeseries(config)

if config.getboolean('sst_timeseries','generate'):
    print ""
    print "Plotting SST time series..."
    sst_timeseries(config)

if config.getboolean('nino34_timeseries','generate'):
    print ""
    print "Plotting Nino3.4 time series..."
    #nino34_timeseries(config)

if config.getboolean('mht_timeseries','generate'):
    print ""
    print "Plotting Meridional Heat Transport (MHT)..."
    #mht_timeseries(config)


if config.getboolean('moc_timeseries','generate'):
    print ""
    print "Plotting Meridional Overturning Circulation (MOC)..."
    #moc_timeseries(config)
    
if config.getboolean('sst_modelvsobs','generate'):
    print ""
    print "Plotting 2-d maps of SST climatologies..."
    sst_modelvsobs(config)


#GENERATE SEA-ICE DIAGNOSTICS
if config.getboolean('seaice_timeseries','generate'):
    print ""
    print "Plotting sea-ice area and volume time series..."
    seaice_timeseries(config)

if config.getboolean('seaice_modelvsobs','generate'):
    print ""
    print "Plotting 2-d maps of sea-ice concentration and thickness climatologies..."
    seaice_modelvsobs(config)


#GENERATE LAND-ICE DIAGNOSTICS


if config.getboolean('plot','displayToScreen'):
   plt.show()

