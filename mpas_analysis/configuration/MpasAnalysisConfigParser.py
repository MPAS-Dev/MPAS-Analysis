import numbers
import ast

"""
A configuratin parser class for MPAS analysis.  MpasAnalysisConfigParser adds
the capabilities to get an option including a default value
(`getWithDefault(section, option, default, ...)`) and to get options
that are lists, tuples, dicts, etc (`getExpression(section, option)`).

Author: Xylar Asay-Davis
Last Modified: 12/03/2016
"""

from ConfigParser import ConfigParser


class MpasAnalysisConfigParser(ConfigParser):

    def getWithDefault(self, section, option, default):
        """
        Get an option, using the supplied default value if the option is not
        present.  The type of `default` is used to determine the type of the
        retrun value, regardless of whether `default` is returned. If `default`
        is a list, tuple, or dict, `getExpression(...)` is used if the option
        is present in the config file.

        Author: Xylar Asay-Davis
        Last Modified: 12/03/2016
        """
        if self.has_section(section):
            if self.has_option(section, option):
                if isinstance(default, numbers.Integral):
                    return self.getint(section, option)
                elif isinstance(default, numbers.Real):
                    return self.getfloat(section, option)
                elif isinstance(default, bool):
                    return self.getboolean(section, option)
                elif isinstance(default, (list, tuple, dict)):
                    return self.getExpression(section, option)
                else:
                    return self.get(section, option)

        # we didn't find the entry so set it to the default
        self.set(section, option, str(default))
        return default

    def getExpression(self, section, option):
        """
        Get an option as an expression (typically a list, though tuples and
        dicts should also work).  `section` and `option` work as in `get(...)`.
        The expression is required to have valid python syntax, so that
        string entries are required to be in single or double quotes.

        Author: Xylar Asay-Davis
        Last Modified: 12/03/2016
        """
        expressionString = self.get(section, option)
        return ast.literal_eval(expressionString)
