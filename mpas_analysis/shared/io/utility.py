"""
IO utility functions

Phillip J. Wolfram, Xylar Asay-Davis
"""

import glob
import os
import random
import string


def paths(*args): # {{{
    """
    Returns glob'd paths in list for arbitrary number of function arguments.
    Note, each expanded set of paths is sorted.

    Phillip J. Wolfram
    """
    paths = []
    for aargs in args:
        paths += sorted(glob.glob(aargs))
    return paths # }}}


def fingerprint_generator(size=12,
                          chars=string.ascii_uppercase + string.digits): # {{{
    """
    Returns a random string that can be used as a unique fingerprint

    Reference:
    http://stackoverflow.com/questions/2257441/random-string-generation-with-upper-case-letters-and-digits-in-python

    Phillip J. Wolfram
    """
    return ''.join(random.choice(chars) for _ in range(size)) # }}}


def make_directories(path):  # {{{
    """
    Make the given path if it does not already exist.

    Returns the path unchanged.

    Author: Xylar Asay-Davis
    """

    try:
        os.makedirs(path)
    except OSError:
        pass
    return path  # }}}


def build_config_full_path(config, section, relativePathOption,
                           relativePathSection=None,
                           defaultPath=None): # {{{
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
    return fullPath # }}}


def check_path_exists(path):  # {{{
    """
    Raise an exception if the given path does not exist.

    Parameters
    ----------
    path : str
        Absolute path

    Raises
    ------
    OSError
        If the path does not exist

    Authors
    -------
    Xylar Asay-Davis
    """
    if not (os.path.isdir(path) or os.path.isfile(path)):
        raise OSError('Path {} not found'.format(path))  # }}}


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
