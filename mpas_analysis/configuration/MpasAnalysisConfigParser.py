import numbers

from ConfigParser import ConfigParser

class MpasAnalysisConfigParser(ConfigParser):

    def getWithDefault(self, section, option, default, listType=None):
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
            return self.get(section, option, str(default))
        return default

    def getlist(self, section, option, listType=str):
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
        listString = ', '.join(str(entry) for entry in value)
        listString = '[%s]'%listString
        self.set(section, option, listString)


