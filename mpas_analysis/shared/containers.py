#!/usr/bin/env python
# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2018 Los Alamos National Security, LLC. All rights reserved.
# Copyright (c) 2018 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2018 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
    Module of custom container data types

    Phillip J. Wolfram
    10/22/2016
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import collections


class ReadOnlyDict(collections.Mapping):  # {{{
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
# }}}

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
