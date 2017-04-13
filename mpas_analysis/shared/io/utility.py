"""
IO utility functions

Phillip J. Wolfram, Xylar Asay-Davis

Last Modified: 03/23/2017
"""

import glob
import os


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


def make_directories(path):  # {{{
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


def build_config_full_path(config, section, relativePathOption,
                           relativePathSection=None,
                           defaultPath=None):
    """
    Returns a full path from a base directory and a relative path

    Parameters
    ----------
    config : MpasAnalysisConfigParser object
        configuration from which to read the path

    section : str
        the name of a section in `config`, which must have an option
        `baseDirectory`

    relativePathOption : str
        the name of an option in `section` of the relative path within
        `baseDirectory` (or possibly an absolute path)

    relativePathSection : str, optional
        the name of a section for `relativePathOption` if not `section`

    defaultPath : str, optional
        the name of a path to return if the resulting path doesn't exist.

    Authors
    -------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/23/2017
    """
    if relativePathSection is None:
        relativePathSection = section

    subDirectory = config.get(relativePathSection, relativePathOption)
    if os.path.isabs(subDirectory):
        fullPath = subDirectory
    else:
        fullPath = '{}/{}'.format(config.get(section, 'baseDirectory'),
                                  subDirectory)

    if defaultPath is not None and not os.path.exists(fullPath):
        fullPath = defaultPath
    return fullPath

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
