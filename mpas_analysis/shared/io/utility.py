#!/usr/bin/env python
"""
IO utility functions

Phillip J. Wolfram, Xylar Asay-Davis

Last Modified: 02/02/2017
"""

import glob


def paths(*args):
    """
    Returns glob'd paths in list for arbitrary number of function arguments.
    Note, each expanded set of paths is sorted.

    Phillip J. Wolfram
    10/25/2016
    """
    paths = []
    for aargs in args:
        paths += sorted(glob.glob(aargs))
    return paths


def buildConfigFullPath(config, section, relativePathOption):
    """
    Returns a full path from a base directory and a relative path

    `config` is an instance of a ConfigParser

    `section` is the name of a section in `config`, which must have an
    option `baseDirectory`

    `relativePathOption` is the name of an option in `section` that
    points to the name of a relative path within `baseDirectory`

    Xylar Asay-Davis

    Last Modified: 02/02/2017
    """
    return '{}/{}'.format(config.get(section, 'baseDirectory'),
                          config.get(section, relativePathOption))

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
