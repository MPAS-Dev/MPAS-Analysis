"""
Unit test infrastructure, adapted from approach of xarray.

Phillip J. Wolfram
10/07/2016
"""

import warnings
from contextlib import contextmanager

import os
from distutils import dir_util
from pytest import fixture

try:
    import unittest2 as unittest
except ImportError:
    import unittest

try:
    import lxml
    has_lxml = True
except ImportError:
    has_lxml = False

try:
    import numpy as np
    has_numpy = True
except ImportError:
    has_numpy = False


def requires_lxml(test):
    return test if has_lxml else unittest.skip('requires lxml')(test)


def requires_numpy(test):
    return test if has_numpy else unittest.skip('requires numpy')(test)


# Adapted from
# http://stackoverflow.com/questions/29627341/pytest-where-to-store-expected-data
@fixture
def loaddatadir(request, tmpdir):
    '''
    Fixture responsible for searching a folder with the same name of test
    module and, if available, moving all contents to a temporary directory so
    tests can use them freely.
    '''
    filename = request.module.__file__
    test_dir, _ = os.path.splitext(filename)

    if os.path.isdir(test_dir):
        dir_util.copy_tree(test_dir, bytes(tmpdir))

    request.cls.datadir = tmpdir


class TestCase(unittest.TestCase):
    def assertEqual(self, a1, a2):
        assert a1 == a2 or (a1 != a1 and a2 != a2)

    def assertLessThan(self, a1, a2):
        assert a1 <= a2

    def assertGreaterThan(self, a1, a2):
        assert a1 >= a2

    @requires_numpy
    def assertArrayEqual(self, a1, a2):
        np.testing.assert_array_equal(a1, a2)

    @requires_numpy
    def assertApproxEqual(self, a1, a2, rtol=1e-5, atol=1e-8):
        assert np.isclose(a1, a2, rtol=rtol, atol=atol)

    @requires_numpy
    def assertArrayApproxEqual(self, a1, a2, rtol=1e-5, atol=1e-8):
        assert np.all(np.isclose(a1, a2, rtol=rtol, atol=atol))

    @contextmanager
    def assertWarns(self, message):
        with warnings.catch_warnings(record=True) as w:
            warnings.filterwarnings('always', message)
            yield
            assert len(w) > 0
            assert all(message in str(wi.message) for wi in w)


