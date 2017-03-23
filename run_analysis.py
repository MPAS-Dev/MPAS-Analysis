#!/usr/bin/env python

"""
Runs MPAS-Analysis via a configuration file (e.g. `config.analysis`)
specifying analysis options.

Author: Xylar Asay-Davis, Phillip J. Wolfram
Last Modified: 03/23/2017
"""

import os
import matplotlib as mpl
import argparse

from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


def checkPathExists(path):  # {{{
    """
    Raise an exception if the given path does not exist.

    Author: Xylar Asay-Davis
    Last Modified: 02/02/2017
    """
    if not (os.path.isdir(path) or os.path.isfile(path)):
        raise OSError('Path {} not found'.format(path))
# }}}


def makeDirectories(path):  # {{{
    """
    Make the given path if it does not already exist.

    Returns the path unchanged.

    Author: Xylar Asay-Davis
    Last Modified: 02/02/2017
    """

    try:
        os.makedirs(path)
    except OSError:
        pass
    return path  # }}}


def checkGenerate(config, analysisName, mpasCore, analysisCategory=None):
    # {{{
    """
    determine if a particular analysis of a particular core and (optionally)
    category should be generated.

    Author: Xylar Asay-Davis
    Last Modified: 02/02/2017
    """
    generateList = config.getExpression('output', 'generate')
    generate = False
    for element in generateList:
        if '_' in element:
            (prefix, suffix) = element.split('_', 1)
        else:
            prefix = element
            suffix = None

        if prefix == 'all':
            if (suffix in [mpasCore, analysisCategory]) or (suffix is None):
                generate = True
        elif prefix == 'no':
            if suffix in [analysisName, mpasCore, analysisCategory]:
                generate = False
        elif element == analysisName:
            generate = True

    return generate  # }}}


def analysis(config):  # {{{
    # set default values of start and end dates for climotologies and
    # timeseries and indices
    for section in ['climatology', 'timeSeries', 'index']:
        startDate = '{:04d}-01-01_00:00:00'.format(
            config.getint(section, 'startYear'))
        if not config.has_option(section, 'startDate'):
            config.set(section, 'startDate', startDate)
        endDate = '{:04d}-12-31_23:59:59'.format(
            config.getint(section, 'endYear'))
        if not config.has_option(section, 'endDate'):
            config.set(section, 'endDate', endDate)

    # Checks on directory/files existence:
    if config.get('runs', 'preprocessedReferenceRunName') != 'None':
        checkPathExists(config.get('oceanPreprocessedReference',
                                   'baseDirectory'))
        checkPathExists(config.get('seaIcePreprocessedReference',
                                   'baseDirectory'))

    generateTimeSeriesSeaIce = checkGenerate(
        config, analysisName='timeSeriesSeaIceAreaVol',  mpasCore='seaIce',
        analysisCategory='timeSeries')
    compareTimeSeriesSeaIceWithObservations = config.getboolean(
            'timeSeriesSeaIceAreaVol', 'compareWithObservations')
    generateRegriddedSeaIce = checkGenerate(
        config, analysisName='regriddedSeaIceConcThick',  mpasCore='seaIce',
        analysisCategory='regriddedHorizontal')

    if ((generateTimeSeriesSeaIce and
         compareTimeSeriesSeaIceWithObservations) or generateRegriddedSeaIce):
        # we will need sea-ice observations.  Make sure they're there
        baseDirectory = config.get('seaIceObservations', 'baseDirectory')
        for observationName in ['areaNH', 'areaSH', 'volNH', 'volSH']:
            fileName = config.get('seaIceObservations', observationName)
            if fileName.lower() == 'none':
                continue
            checkPathExists('{}/{}'.format(baseDirectory, fileName))

    # choose the right rendering backend, depending on whether we're displaying
    # to the screen
    if not config.getboolean('plot', 'displayToScreen'):
        mpl.use('Agg')
    import matplotlib.pyplot as plt

    # analysis can only be imported after the right MPL renderer is selected

    # GENERATE OCEAN DIAGNOSTICS
    if checkGenerate(config, analysisName='timeSeriesOHC', mpasCore='ocean',
                     analysisCategory='timeSeries'):
        print ""
        print "Plotting OHC time series..."
        from mpas_analysis.ocean.ohc_timeseries import ohc_timeseries
        ohc_timeseries(config)

    if checkGenerate(config, analysisName='timeSeriesSST', mpasCore='ocean',
                     analysisCategory='timeSeries'):
        print ""
        print "Plotting SST time series..."
        from mpas_analysis.ocean.sst_timeseries import sst_timeseries
        sst_timeseries(config)

    if checkGenerate(config, analysisName='indexNino34',
                     mpasCore='ocean', analysisCategory='index'):
        print ""
        print "Plotting Nino3.4 time series and power spectrum...."
        from mpas_analysis.ocean.nino34_index import nino34_index
        nino34_index(config)

#    if checkGenerate(config, analysisName='timeSeriesMHT', mpasCore='ocean',
#                     analysisCategory='timeSeries'):
#        print ""
#        print "Plotting Meridional Heat Transport (MHT)..."
#        from mpas_analysis.ocean.mht_timeseries import mht_timeseries
#        mht_timeseries(config)

    if checkGenerate(config, analysisName='regriddedSST', mpasCore='ocean',
                     analysisCategory='regriddedHorizontal'):
        print ""
        print "Plotting 2-d maps of SST climatologies..."
        from mpas_analysis.ocean.ocean_modelvsobs import ocn_modelvsobs
        ocn_modelvsobs(config, 'sst')

    if checkGenerate(config, analysisName='regriddedMLD', mpasCore='ocean',
                     analysisCategory='regriddedHorizontal'):
        print ""
        print "Plotting 2-d maps of MLD climatologies..."
        from mpas_analysis.ocean.ocean_modelvsobs import ocn_modelvsobs
        ocn_modelvsobs(config, 'mld')

    if checkGenerate(config, analysisName='regriddedSSS', mpasCore='ocean',
                     analysisCategory='regriddedHorizontal'):
        print ""
        print "Plotting 2-d maps of SSS climatologies..."
        from mpas_analysis.ocean.ocean_modelvsobs import ocn_modelvsobs
        ocn_modelvsobs(config, 'sss')

    if checkGenerate(config, analysisName='streamfunctionMOC',
                     mpasCore='ocean',
                     analysisCategory='streamfunctionMOC'):
        print ""
        print "Plotting streamfunction of Meridional Overturning Circulation (MOC)..."
        from mpas_analysis.ocean.meridional_overturning_circulation \
            import moc_streamfunction
        moc_streamfunction(config)

    # GENERATE SEA-ICE DIAGNOSTICS
    if checkGenerate(config, analysisName='timeSeriesSeaIceAreaVol',
                     mpasCore='seaIce', analysisCategory='timeSeries'):
        print ""
        print "Plotting sea-ice area and volume time series..."
        from mpas_analysis.sea_ice.timeseries import seaice_timeseries
        seaice_timeseries(config)

    if checkGenerate(config, analysisName='regriddedSeaIceConcThick',
                     mpasCore='seaIce',
                     analysisCategory='regriddedHorizontal'):
        print ""
        print "Plotting 2-d maps of sea-ice concentration and thickness " \
            "climatologies..."
        from mpas_analysis.sea_ice.modelvsobs import seaice_modelvsobs
        seaice_modelvsobs(config)

    # GENERATE LAND-ICE DIAGNOSTICS

    if config.getboolean('plot', 'displayToScreen'):
        plt.show()

    return  # }}}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-g", "--generate", dest="generate",
                        help="A list of analysis modules to generate "
                        "(nearly identical generate option in config file).",
                        metavar="ANALYSIS1[,ANALYSIS2,ANALYSIS3,...]")
    parser.add_argument('configFiles', metavar='CONFIG',
                        type=str, nargs='+', help='config file')
    args = parser.parse_args()

    # add config.default to cover default not included in the config files
    # provided on the command line
    defaultConfig = '{}/config.default'.format(
        os.path.dirname(os.path.realpath(__file__)))
    if os.path.exists(defaultConfig):
        configFiles = [defaultConfig] + args.configFiles
    else:
        print 'WARNING: Did not find config.default.  Assuming other config ' \
            'file(s) contain a\nfull set of configuration options.'
        configFiles = args.configFiles

    config = MpasAnalysisConfigParser()
    config.read(configFiles)

    if args.generate:
        # overwrite the 'generate' in config with a string that parses to
        # a list of string
        generateList = args.generate.split(',')
        generateString = ', '.join(["'{}'".format(element)
                                    for element in generateList])
        generateString = '[{}]'.format(generateString)
        config.set('output', 'generate', generateString)

    analysis(config)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
