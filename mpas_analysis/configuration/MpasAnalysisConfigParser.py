import numbers

"""
A configuratin parser class for MPAS analysis.  MpasAnalysisConfigParser adds
the capabilities to get an option including a default value
(`getWithDefault(section, option, default, ...)`) and to get and set options
that are lists (`getList(section, option, ...)` and
`setList(secton, option, value)`).

Author: Xylar Asay-Davis
Last Modified: 12/03/2016
"""

from ConfigParser import ConfigParser


class MpasAnalysisConfigParser(ConfigParser):

    def getWithDefault(self, section, option, default, listType=None):
        """
        Get an option, using the supplied default value if the option is not
        present.  The type of `default` is used to determine the type of the
        retrun value, regardless of whether `default` is returned. If `default`
        is a list and the given option is present, `listType` may be supplied
        to determine the type of the resulting list.  If `listType` is not
        supplied and `default` is a list, the default type of entries in the
        resulting list is the type of the first item in `default`.

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
                elif isinstance(default, list):
                    if (listType is None) and (len(default) > 0):
                        # use the type of default[0] for listType
                        listType = type(default[0])

                    return self.getlist(section, option, listType)
                else:
                    return self.get(section, option)

        # we didn't find the entry so set it to the default
        if isinstance(default, list):
            self.setlist(section, option, default)
        else:
            self.set(section, option, str(default))
        return default

    def getlist(self, section, option, listType=str):
        """
        Get an option as a list.  `section` and `option` work as in `get(...)`.
        `listType` is the type of each entry in the resulting list (`str` by
        default).  Quotes around strings in the list will be stripped (but
        are not required).

        Author: Xylar Asay-Davis
        Last Modified: 12/03/2016
        """
        listString = self.get(section, option)
        # first, strip the brackets
        listString = listString.strip('[]')
        if (listType == str) and ('"' in listString):
            # assume all entries are strings separated by quotes,
            # but every other one is the comma separating list entries
            stringList = listString.split('"')[0::2]
        elif (listType == str) and ("'" in listString):
            # assume all entries are strings separated by quotes,
            # but every other one is the comma separating list entries
            stringList = listString.split("'")[0::2]
        else:
            stringList = listString.split(',')
        outList = []
        for string in stringList:
            outList.append(listType(string.strip()))
        return outList

    def setlist(self, section, option, value):
        """
        Sets an option to be a list based on entries in `value`, which must
        be iterable.  Entries are converted to strings with `str`, so
        it may be desirable to convert floating-point values to stirngs before
        calling this routine if more constrol over the precision is required.

        Author: Xylar Asay-Davis
        Last Modified: 12/03/2016
        """
        listString = ', '.join(str(entry) for entry in value)
        listString = '[{}]'.format(listString)
        self.set(section, option, listString)
