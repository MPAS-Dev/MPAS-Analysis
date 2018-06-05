#!/usr/bin/env python
# Copyright (c) 2017,  Los Alamos National Security, LLC (LANS)
# and the University Corporation for Atmospheric Research (UCAR).
#
# Unless noted otherwise source code is licensed under the BSD license.
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at http://mpas-dev.github.com/license.html
#

"""
A script for downloading the input data set from public repository for
MPAS-Analysis to work. The input data set includes: pre-processed
observations data, MPAS mapping files and MPAS regional mask files
(which are used for the MOC computation), for a subset of MPAS meshes.
"""
# Authors
# -------
# Milena Veneziani
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals


import os
import argparse
import pkg_resources

from mpas_analysis.shared.io.download import download_files


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--outDir", dest="outDir", required=True,
                        help="Directory where MPAS-Analysis input data will"
                             "be downloaded")
    args = parser.parse_args()

    try:
        os.makedirs(args.outDir)
    except OSError:
        pass

    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics'
    analysisFileList = pkg_resources.resource_string(
            'mpas_analysis', 'obs/analysis_input_files').decode('utf-8')

    # remove any empty strings from the list
    analysisFileList = list(filter(None, analysisFileList.split('\n')))
    download_files(analysisFileList, urlBase, args.outDir)
