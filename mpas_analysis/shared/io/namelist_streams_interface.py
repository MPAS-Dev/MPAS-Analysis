#!/usr/bin/env python
"""
Module of classes/routines to manipulate fortran namelist and streams
files.

Authors
-------
Phillip Wolfram, Xylar Asay-Davis

Last modified
-------------
02/06/2017
"""

from lxml import etree
import re
import os.path

from ..containers import ReadOnlyDict
from .utility import paths
from ..timekeeping.utility import string_to_datetime, string_to_relative_delta


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

    Authors
    -------
    Phillip Wolfram, Xylar Asay-Davis

    Last modified
    -------------
    02/06/2017
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
        """ Accessor for bracket noation, e.g., nml['field'], returns string
        """
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
    # }}}


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

    def readpath(self, streamName, startDate=None, endDate=None,
                 calendar=None):
        """
        Given the name of a stream and optionally start and end dates and a
        calendar type, returns a list of files that match the file template in
        the stream.

        Parameters
        ----------
        streamName : string
            The name of a stream that produced the files

        startDate, endDate : string or datetime.datetime, optional
            String or datetime.datetime objects identifying the beginning
            and end dates to be found.

            Note: a buffer of one output interval is subtracted from startDate
            and added to endDate because the file date might be the first
            or last date contained in the file (or anything in between).

        calendar: {'gregorian', 'gregorian_noleap'}, optional
            The name of one of the calendars supported by MPAS cores, and is
            required if startDate and/or endDate are supplied

        Returns
        -------
        fileList : list
            A list of file names produced by the stream that fall between
            the startDate and endDate (if supplied)

        Raises
        ------
        ValueError
            If no files from the stream are found.

        Author
        ------
        Xylar Asay-Davis

        Last modified
        -------------
        02/04/2017
        """
        template = self.read(streamName, 'filename_template')
        if template is None:
            raise ValueError('Stream {} not found in streams file {}.'.format(
                streamName, self.fname))
        replacements = {'$Y': '[0-9][0-9][0-9][0-9]',
                        '$M': '[0-9][0-9]',
                        '$D': '[0-9][0-9]',
                        '$S': '[0-9][0-9][0-9][0-9][0-9]',
                        '$h': '[0-9][0-9]',
                        '$m': '[0-9][0-9]',
                        '$s': '[0-9][0-9]'}

        path = template
        for old in replacements:
            path = path.replace(old, replacements[old])

        if not os.path.isabs(path):
            # this is not an absolute path, so make it an absolute path
            path = '{}/{}'.format(self.streamsdir, path)

        fileList = paths(path)

        if len(fileList) == 0:
            raise ValueError(
                "Path {} in streams file {} for '{}' not found.".format(
                    path, self.fname, streamName))

        if (startDate is None) and (endDate is None):
            return fileList

        output_interval = self.read(streamName, 'output_interval')
        if output_interval is None:
            # There's no file interval, so hard to know what to do
            # let's put a buffer of a year on each side to be safe
            offsetDate = string_to_relative_delta(dateString='0001-00-00',
                                                  calendar=calendar)
        else:
            offsetDate = string_to_relative_delta(dateString=output_interval,
                                                  calendar=calendar)

        if startDate is not None:
            # read one extra file before the start date to be on the safe side
            if isinstance(startDate, str):
                startDate = string_to_datetime(startDate)
            try:
                startDate -= offsetDate
            except (ValueError, OverflowError):
                # if the startDate would be out of range after subtracting
                # the offset, we'll stick with the starDate as it is
                pass

        if endDate is not None:
            # read one extra file after the end date to be on the safe side
            if isinstance(endDate, str):
                endDate = string_to_datetime(endDate)
            try:
                endDate += offsetDate
            except (ValueError, OverflowError):
                # if the endDate would be out of range after adding
                # the offset, we'll stick with the endDate as it is
                pass

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
            fileDate = string_to_datetime(fileDateString)
            add = True
            if startDate is not None and startDate > fileDate:
                add = False
            if endDate is not None and endDate < fileDate:
                add = False
            if add:
                outFileList.append(fileName)

        return outFileList

    def has_stream(self, streamName):
        """
        Returns True if the streams file has a stream with the given
        streamName, otherwise returns False.

        Xylar Asay-Davis
        Last modified: 12/04/2016
        """
        for stream in self.root:
            # assumes streamname is unique in XML
            if stream.get('name') == streamName:
                return True
        return False

    def find_stream(self, possibleStreams):
        """
        If one (or more) of the names in possibleStreams is a stream in this
        streams file, returns the first match.  If no match is found, raises
        a ValueError.

        Xylar Asay-Davis
        Last modified: 12/07/2016
        """
        for streamName in possibleStreams:
            if self.has_stream(streamName):
                return streamName

        raise ValueError('Stream {} not found in streams file {}.'.format(
            streamName, self.fname))

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
