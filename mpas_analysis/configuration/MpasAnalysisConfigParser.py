"""
A configuratin parser class for MPAS analysis.  MpasAnalysisConfigParser adds
the capabilities to get an option including a default value
(`getWithDefault(section, option, default, ...)`) and to get options
that are lists, tuples, dicts, etc (`getExpression(section, option)`).

Author: Xylar Asay-Davis, Phillip J. Wolfram
Last Modified: 02/27/2017
"""

import numbers
import ast
import numpy as np
from ConfigParser import ConfigParser


npallow = dict(linspace=np.linspace, xrange=xrange, range=range,
               arange=np.arange, pi=np.pi, Pi=np.pi, __builtins__=None)


class MpasAnalysisConfigParser(ConfigParser):

    def getWithDefault(self, section, option, default):
        """
        Get an option, using the supplied default value if the option is not
        present.  The type of `default` is used to determine the type of the
        retrun value, regardless of whether `default` is returned. If `default`
        is a list, tuple, or dict, `getExpression(...)` is used if the option
        is present in the config file.

        Author: Xylar Asay-Davis
        Last Modified: 02/27/2017
        """
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
        dicts should also work).  `section` and `option` work as in `get(...)`.
        The expression is required to have valid python syntax, so that
        string entries are required to be in single or double quotes.

        If `elementType` is supplied, each element in a list or tuple, or each
        value in a dictionary are cast to this type.  This is likely most
        useful for ensuring that all elements in a list of numbers are of type
        float, rather than int, when the distinction is important.

        If `usenumpyfunc` is True, expression is evaluated within the context
        of having selected numpy and / or np functionality available.

        Author: Xylar Asay-Davis, Phillip J. Wolfram
        Last Modified: 01/31/2017
        """
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
