#!/usr/bin/env python

"""
Runs MPAS-Analysis via a configuration file (e.g. `config.analysis`)
specifying analysis options.

Authors
-------
Xylar Asay-Davis, Phillip J. Wolfram
"""

import os
import matplotlib as mpl
import argparse
import traceback
import sys
import warnings

from mpas_analysis.configuration.MpasAnalysisConfigParser \
    import MpasAnalysisConfigParser


def update_generate(config, generate):  # {{{
    """
    Update the 'generate' config option using a string from the command line.

    Author: Xylar Asay-Davis
    Last Modified: 03/07/2017
    """

    # overwrite the 'generate' in config with a string that parses to
    # a list of string
    generateList = generate.split(',')
    generateString = ', '.join(["'{}'".format(element)
                                for element in generateList])
    generateString = '[{}]'.format(generateString)
    config.set('output', 'generate', generateString)  # }}}


def build_analysis_list(config):  # {{{
    """
    Build a list of analysis modules based on the 'generate' config option.

    Author: Xylar Asay-Davis
    Last Modified: 03/07/2017
    """

    # choose the right rendering backend, depending on whether we're displaying
    # to the screen
    if not config.getboolean('plot', 'displayToScreen'):
        mpl.use('Agg')

    # analysis can only be imported after the right MPL renderer is selected
    from mpas_analysis import ocean as ocean_tasks
    from mpas_analysis import sea_ice as sea_ice_tasks

    # analyses will be a list of analysis classes
    analyses = []

    # Ocean Analyses
    analyses.append(ocean_tasks.TimeSeriesOHC(config))
    # analyses.append(ocean_tasks.TimeSeriesSST(config))
    # analyses.append(ocean_tasks.IndexNino34(config))
    # analyses.append(ocean_tasks.MeridionalHeatTransport(config))
    # analyses.append(ocean_tasks.MeridionalOverturningCirculation(config))

    analyses.append(ocean_tasks.ClimatologyMapSST(config))
    analyses.append(ocean_tasks.ClimatologyMapMLD(config))
    analyses.append(ocean_tasks.ClimatologyMapSSS(config))

    # Sea Ice Analyses
    # analyses.append(sea_ice_tasks.TimeSeriesSeaIce(config))
    analyses.append(sea_ice_tasks.ClimatologyMapSeaIce(config))

    # check which analysis we actually want to generate and only keep those
    analysesToGenerate = []
    for analysisTask in analyses:
        # for each anlaysis module, check if we want to generate this task
        # and if the analysis task has a valid configuration
        if analysisTask.check_generate():
            add = False
            try:
                analysisTask.setup_and_check()
                add = True
            except:
                traceback.print_exc(file=sys.stdout)
                print "ERROR: analysis module {} failed during check and " \
                    "will not be run".format(analysisTask.taskName)
            if add:
                analysesToGenerate.append(analysisTask)

    return analysesToGenerate  # }}}


def run_analysis(config, analyses):  # {{{

    # run each analysis task
    for analysisTask in analyses:
        try:
            analysisTask.run()
        except:
            traceback.print_exc(file=sys.stdout)
            print "ERROR: analysis module {} failed during run".format(
                analysisTask.taskName)

    if config.getboolean('plot', 'displayToScreen'):
        import matplotlib.pyplot as plt
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
        warnings.warn('WARNING: Did not find config.default.  Assuming other '
                      'config file(s) contain a\n'
                      'full set of configuration options.')
        configFiles = args.configFiles

    config = MpasAnalysisConfigParser()
    config.read(configFiles)

    if args.generate:
        update_generate(config, args.generate)

    analyses = build_analysis_list(config)

    run_analysis(config, analyses)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
