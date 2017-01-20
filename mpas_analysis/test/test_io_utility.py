"""
Unit test infrastructure, adapted from approach of xarray.

Phillip J. Wolfram
10/07/2016
"""

import os
import pytest
from mpas_analysis.test import TestCase, loaddatadir
from mpas_analysis.shared.io import paths


@pytest.mark.usefixtures("loaddatadir")
class TestPaths(TestCase):
    def test_paths(self):
        os.chdir(bytes(self.datadir))
        self.assertEquals(paths('[0-9]*', '[a-z]*'),
                          ['0.txt', '1.txt', '2.txt', 'a.txt', 'b.txt',
                           'c.txt'])


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
