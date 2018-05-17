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
import requests
import progressbar


# From https://stackoverflow.com/a/1094933/7728169
def sizeof_fmt(num, suffix='B'):
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


def download_analysis_files(outDir):
    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics'
    analysisFileList = pkg_resources.resource_string(
            'mpas_analysis', 'obs/analysis_input_files').decode('utf-8')

    for fileName in analysisFileList.split('\n'):
        outFileName = '{}/{}'.format(outDir, fileName)
        # outFileName contains full path, so we need to make the relevant
        # subdirectories if they do not exist already
        directory = os.path.dirname(outFileName)
        try:
            os.makedirs(directory)
        except OSError:
            pass

        url = '{}/{}'.format(urlBase, fileName)
        try:
            response = requests.get(url, stream=True)
            totalSize = response.headers.get('content-length')
        except requests.exceptions.RequestException:
            print('  {} could not be reached!'.format(fileName))
            continue

        if totalSize is None:
            # no content length header
            if not os.path.exists(outFileName):
                with open(outFileName, 'wb') as f:
                    print('Downloading {}...'.format(fileName))
                    try:
                        f.write(response.content)
                    except requests.exceptions.RequestException:
                        print('  {} failed!'.format(fileName))
                    else:
                        print('  {} done.'.format(fileName))
        else:
            # we can do the download in chunks and use a progress bar, yay!

            totalSize = int(totalSize)
            if os.path.exists(outFileName) and \
                    totalSize == os.path.getsize(outFileName):
                # we already have the file, so just continue
                continue

            print('Downloading {} ({})...'.format(fileName,
                                                  sizeof_fmt(totalSize)))
            widgets = [progressbar.Percentage(), ' ', progressbar.Bar(),
                       ' ', progressbar.ETA()]
            bar = progressbar.ProgressBar(widgets=widgets,
                                          maxval=totalSize).start()
            size = 0
            with open(outFileName, 'wb') as f:
                try:
                    for data in response.iter_content(chunk_size=4096):
                        size += len(data)
                        f.write(data)
                        bar.update(size)
                    bar.finish()
                except requests.exceptions.RequestException:
                    print('  {} failed!'.format(fileName))
                else:
                    print('  {} done.'.format(fileName))


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

    download_analysis_files(args.outDir)
