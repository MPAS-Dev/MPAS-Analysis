#!/usr/bin/env python
"""
    Module of custom container data types

    Phillip J. Wolfram
    10/22/2016
"""

import collections
class ReadOnlyDict(collections.Mapping): #{{{
    """ Read only-dictionary
    http://stackoverflow.com/questions/19022868/how-to-make-dictionary-read-only-in-python
    310/22/2016
    """

    def __init__(self, data):
        self._data = data

    # overloads [] operator
    def __getitem__(self, key):
        return self._data[key]

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)
#}}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
