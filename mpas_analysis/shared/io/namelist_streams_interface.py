#!/usr/bin/env python
"""
    Module of classes / routines to manipulate fortran namelist and streams
    files.

    Phillip Wolfram, Xylar Asay-Davis
    Last modified: 11/02/2016
"""

from lxml import etree
import re
import os.path

from ..containers import ReadOnlyDict
from .utility import paths
from ..timekeeping.Date import Date

def convert_namelist_to_dict(fname, readonly=True):
    """
    Converts a namelist file to key-value pairs in dictionary.

    Phillip J Wolfram
    10/22/2016
    """
    # form dictionary
    nml = dict()

    regex = re.compile(r"^\s*(.*?)\s*=\s*['\"]*(.*?)['\"]*\s*\n")
    with open(fname) as f:
        for line in f:
            match = regex.findall(line)
            if len(match) > 0:
                # assumes that there is only one match per line
                nml[match[0][0]] = match[0][1]
    if readonly:
        nml = ReadOnlyDict(nml)

    return nml


class NameList:
    """
    Class for fortran manipulation of namelist files, provides
    read and write functionality

    Phillip Wolfram, Xylar Asay-Davis
    Last modified: 11/02/2016
    """

    # constructor
    def __init__(self, fname, path=None):
        """
        parse the namelist file given by fname. If the optional argument
        path is provided, and fname contains a relative path, fname is
        relative to path, rather than the current working directory.
        """
        if not os.path.isabs(fname) and path is not None:
            # only the file name was given, not the absolute path, and
            # a path was provided, so we will assume the namelist
            # file is actually in the path
            fname = '{}/{}'.format(path, fname)

        # input file name
        self.fname = fname
        # get values
        self.nml = convert_namelist_to_dict(fname)

    # note following accessors do not do type casting
    def __getattr__(self, key):
        """ Accessor for dot noation, e.g., nml.field, returns string """
        return self.nml[key]

    # provide accessor for dictionary notation (returns string)
    def __getitem__(self, key):
        """ Accessor for bracket noation, e.g., nml['field'], returns string """
        return self.nml[key]

    # provide accessors for get, getint, getfloat, getbool with appropriate
    # casting for comparable behavior with config files #{{{
    def get(self, key):
        return self.nml[key]

    def getint(self, key):
        return int(self.nml[key])

    def getfloat(self, key):
        return float(self.nml[key])

    def getbool(self, key):
        if 'True' in self.nml[key] or 'true' in self.nml[key]:
            return True
        else:
            return False
    #}}}


class StreamsFile:
    """
    Class to read in streams configuration file, provdies
    read and write functionality

    Phillip Wolfram, Xylar Asay-Davis
    Last modified: 11/02/2016
    """

    def __init__(self, fname, streamsdir=None):
        """
        parse the streams file given by fname. If the optional argument
        streamsdir is provided, it is the base path to both the output streams
        data and the sreams file (the latter only if fname is a relative path).
        """
        if not os.path.isabs(fname) and streamsdir is not None:
            # only the file name was given, not the absolute path, and
            # a streamsdir was provided, so we will assume the streams
            # file is actually in the streamsdir
            fname = '{}/{}'.format(streamsdir, fname)

        self.fname = fname
        self.xmlfile = etree.parse(fname)
        self.root = self.xmlfile.getroot()

        if streamsdir is None:
            # get the absolute path to the directory where the
            # streams file resides (used to determine absolute paths
            # to file names referred to in streams)
            self.streamsdir = os.path.dirname(os.path.abspath(fname))
        else:
            self.streamsdir = streamsdir

    def read(self, streamname, attribname):
        """ name is a list of name entries terminanting in some value
        """
        for stream in self.root:
            # assumes streamname is unique in XML
            if stream.get('name') == streamname:
                return stream.get(attribname)
        return None

    def readpath(self, streamName, startDate=None, endDate=None):
        """
        Returns a list of files that match the file template in the
        stream streamName with attribute attribName.  If the startDate
        and/or endDate are supplied, only files on or after the starDate and/or
        on or before the endDate are included in the file list.
        """
        template = self.read(streamName, 'filename_template')
        replacements = {'$Y':'[0-9][0-9][0-9][0-9]',
                        '$M':'[0-9][0-9]',
                        '$D':'[0-9][0-9]',
                        '$S':'[0-9][0-9][0-9][0-9][0-9]',
                        '$h':'[0-9][0-9]',
                        '$m':'[0-9][0-9]',
                        '$s':'[0-9][0-9]'}

        path = template
        for old in replacements:
            path = path.replace(old,replacements[old])

        if not os.path.isabs(path):
            # this is not an absolute path, so make it an absolute path
            path = '{}/{}'.format(self.streamsdir, path)

        fileList = paths(path)

        if (startDate is None) and (endDate is None):
            return fileList

        output_interval = self.read(streamName, 'output_interval')
        if output_interval is None:
            # There's no file interval, so hard to know what to do
            # let's put a buffer of a year on each side to be safe
            offsetDate = Date(dateString='0001-00-00', isInterval=True)
        else:
            offsetDate = Date(dateString=output_interval, isInterval=True)

        if startDate is not None:
            # read one extra file before the start date to be on the safe side
            startDate = Date(startDate) - offsetDate
        if endDate is not None:
            # read one extra file after the end date to be on the safe side
            endDate = Date(endDate) + offsetDate

        # remove any path that's part of the template
        template = os.path.basename(template)
        dateStartIndex = template.find('$')
        if dateStartIndex == -1:
            # there is no date in the template, so we can't exclude any files
            # based on date
            return fileList
        dateEndOffset = len(template) - (template.rfind('$')+2)

        outFileList = []
        for fileName in fileList:
            # get just the
            baseName = os.path.basename(fileName)
            dateEndIndex = len(baseName) - dateEndOffset
            fileDateString = baseName[dateStartIndex:dateEndIndex]
            fileDate = Date(fileDateString)
            add = True
            if startDate is not None and startDate > fileDate:
                add = False
            if endDate is not None and endDate < fileDate:
                add = False
            if add:
                outFileList.append(fileName)

        return outFileList

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
