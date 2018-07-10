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
Unit test infrastructure, adapted from approach of xarray.

Phillip J. Wolfram
10/07/2016
"""

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import os
import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.io import paths


@pytest.mark.usefixtures("loaddatadir")
class TestPaths(TestCase):
    def test_paths(self):
        os.chdir(str(self.datadir))
        self.assertEquals(paths('[0-9]*', '[a-z]*'),
                          ['0.txt', '1.txt', '2.txt', 'a.txt', 'b.txt',
                           'c.txt'])


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
