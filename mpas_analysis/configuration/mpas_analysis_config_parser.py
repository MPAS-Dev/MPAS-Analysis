# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2020 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2020 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2020 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE
"""
A configuratin parser class for MPAS analysis.  MpasAnalysisConfigParser adds
the capabilities to get an option including a default value
(``getWithDefault(section, option, default, ...)``) and to get options
that are lists, tuples, dicts, etc (``getExpression(section, option)``).
"""
# Authors
# -------
# Xylar Asay-Davis, Phillip J. Wolfram

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numbers
import ast
import six
import numpy as np

from six.moves.configparser import ConfigParser

if six.PY3:
    xrange = range

npallow = dict(linspace=np.linspace, xrange=xrange, range=range, array=np.array,
               arange=np.arange, pi=np.pi, Pi=np.pi, int=int, __builtins__=None)


class MpasAnalysisConfigParser(ConfigParser):

    def getWithDefault(self, section, option, default):
        """
        Get an option, using the supplied default value if the option is not
        present.

        Parameters
        ----------
        section : str
            The section in the config file

        option : str
            The option in the config file

        default : one of bool, int, float, list, tuple, dict, str
            The default value if the option and/or section is not found, used
            to determine the type of the option if it *is* found. If
            ``default`` is a list, tuple, or dict, ``getExpression(...)`` is
            used if the option is present in the config file.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        if self.has_section(section):
            if self.has_option(section, option):
                if isinstance(default, bool):
                    return self.getboolean(section, option)
                elif isinstance(default, numbers.Integral):
                    return self.getint(section, option)
                elif isinstance(default, numbers.Real):
                    return self.getfloat(section, option)
                elif isinstance(default, (list, tuple, dict)):
                    return self.getExpression(section, option)
                else:
                    return self.get(section, option)

        # we didn't find the entry so set it to the default
        self.set(section, option, str(default))
        return default

    def getExpression(self, section, option, elementType=None,
                      usenumpyfunc=False):
        """
        Get an option as an expression (typically a list, though tuples and
        dicts are also availabe).  The expression is required to have valid
        python syntax, so that string entries are required to be in single or
        double quotes.

        Parameters
        ----------
        section : str
            The section in the config file

        option : str
            The option in the config file

        elementType : (bool, int, float, etc.), optional
            If supplied, each element in a list or tuple, or
            each value in a dictionary are cast to this type.  This is likely
            most useful for ensuring that all elements in a list of numbers are
            of type float, rather than int, when the distinction is important.

        usenumpyfunc : bool, optional
            If ``True``, the expression is evaluated including functionality
            from the numpy package (which can be referenced either as ``numpy``
            or ``np``).
        """
        # Authors
        # -------
        # Xylar Asay-Davis, Phillip J. Wolfram

        expressionString = self.get(section, option)
        if usenumpyfunc:
            assert '__' not in expressionString, \
                "'__' is not allowed in {} "\
                "for `usenumpyfunc=True`".format(expressionString)
            sanitizedstr = expressionString.replace('np.', '')\
                                           .replace('numpy.', '')\
                                           .replace('__', '')
            result = eval(sanitizedstr, npallow)
        else:
            result = ast.literal_eval(expressionString)

        if elementType is not None:
            if isinstance(result, (list, tuple)):
                result = [elementType(element) for element in result]
            elif isinstance(result, dict):
                for key in result:
                    result[key] = elementType(result[key])

        return result
