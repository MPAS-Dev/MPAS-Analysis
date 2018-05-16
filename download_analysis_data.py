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

from __future__ import absolute_import, division, print_function, \
    unicode_literals


import os
import threading
import argparse
import pkg_resources
try:
    # python 3
    from queue import Queue
except ImportError:
    # python 2
    from Queue import Queue


try:
    # python 3
    from urllib.request import urlretrieve
    from urllib.error import HTTPError
except ImportError:
    # python 2
    from urllib import urlretrieve
    from urllib2 import HTTPError


def download_analysis_files(outDir):
    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics'
    analysisFileList = pkg_resources.resource_string(
            'mpas_analysis', 'obs/analysis_input_files').decode('utf-8')

    fileList = []
    for fileName in analysisFileList.split('\n'):
        fileList.append((urlBase, fileName))

    queue = Queue()
    # Set up some threads to fetch the data
    threadCount = 2
    for index in range(threadCount):
        worker = threading.Thread(target=download_worker,
                                  args=(queue, outDir))
        worker.setDaemon(True)
        worker.start()

    for fileData in fileList:
        queue.put(fileData)
    queue.join()


def download_worker(queue, outDir):
    while True:
        urlBase, fileName = queue.get()
        outFileName = '{}/{}'.format(outDir, fileName)
        # outFileName contains full path, so we need to make the relevant
        # subdirectories if they do not exist already
        directory = os.path.dirname(outFileName)
        try:
            os.makedirs(directory)
        except OSError:
            pass

        if not os.path.exists(outFileName):
            print('Downloading {}...'.format(fileName))
            try:
                urlretrieve('{}/{}'.format(urlBase, fileName), outFileName)
            except HTTPError:
                print('  {} failed!'.format(fileName))
            else:
                print('  {} done.'.format(fileName))
        queue.task_done()


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
