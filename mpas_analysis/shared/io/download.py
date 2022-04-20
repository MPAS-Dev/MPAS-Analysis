#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

"""
Utilities for downloading files
"""
# Authors
# -------
# Milena Veneziani
# Xylar Asay-Davis


import os
import requests
import progressbar


# From https://stackoverflow.com/a/1094933/7728169
def sizeof_fmt(num, suffix='B'):
    """
    Covert a number of bytes to a human-readable file size
    """
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


def download_files(fileList, urlBase, outDir, verify=True):
    """
    Download a list of files from a URL to a directory
    """
    # Authors
    # -------
    # Milena Veneziani
    # Xylar Asay-Davis

    session = requests.Session()
    if not verify:
        session.verify = False

    for fileName in fileList:
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
            response = session.get(url, stream=True)
            totalSize = response.headers.get('content-length')
        except requests.exceptions.RequestException:
            print('  {} could not be reached!'.format(url))
            continue

        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as e:
            print('ERROR while downloading {}:'.format(fileName))
            print(e)
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
                                          max_value=totalSize).start()
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
