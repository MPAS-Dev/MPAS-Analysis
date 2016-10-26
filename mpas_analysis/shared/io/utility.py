#!/usr/bin/env python
"""
IO utility functions

Phillip J. Wolfram
10/25/2016
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

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
