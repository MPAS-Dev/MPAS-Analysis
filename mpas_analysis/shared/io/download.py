#!/usr/bin/env python

"""
Utilities for downloading files
"""
# Authors
# -------
# Milena Veneziani
# Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals


import os
import requests
import progressbar


# From https://stackoverflow.com/a/1094933/7728169
def sizeof_fmt(num, suffix='B'):
    '''
    Covert a number of bytes to a human-readable file size
    '''
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


def download_files(fileList, urlBase, outDir):
    '''
    Download a list of files from a URL to a directory
    '''
    # Authors
    # -------
    # Milena Veneziani
    # Xylar Asay-Davis

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
