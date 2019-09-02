#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

"""
Entry points for downloading data (either for cartopy or for MPAS-Analysis
itself)
"""
# Authors
# -------
# Xylar Asay-Davis, Phillip J. Wolfram, Milena Veneziani

from __future__ import absolute_import, division, print_function, \
    unicode_literals


import argparse
import pkg_resources
import os
import cartopy.io.shapereader as shpreader

from mpas_analysis.shared.io.download import download_files


def download_analysis_data():
    """
    Entry point for downloading the input data set from public repository for
    MPAS-Analysis to work. The input data set includes: pre-processed
    observations data, MPAS mapping files and MPAS regional mask files
    (which are used for the MOC computation), for a subset of MPAS meshes.
    """

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where MPAS-Analysis input data will"
                             "be downloaded")
    parser.add_argument("-d", "--dataset", dest="dataset", default='analysis',
                        help="Directory where MPAS-Analysis input data will"
                             "be downloaded")
    args = parser.parse_args()

    try:
        os.makedirs(args.outDir)
    except OSError:
        pass

    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics'
    analysisFileList = pkg_resources.resource_string(
        'mpas_analysis',
        'obs/{}_input_files'.format(args.dataset)).decode('utf-8')

    # remove any empty strings from the list
    analysisFileList = list(filter(None, analysisFileList.split('\n')))
    download_files(analysisFileList, urlBase, args.outDir, verify=True)


def download_natural_earth_110m():
    for name in ['ocean', 'coastline', 'land']:
        shpfilename = shpreader.natural_earth(resolution='110m',
                                              category='physical',
                                              name=name)
        shpreader.Reader(shpfilename)

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
