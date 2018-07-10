#!/usr/bin/env python
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
Module of classes/routines to manipulate fortran namelist and streams
files.
"""
# Authors
# -------
# Phillip Wolfram, Xylar Asay-Davis

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import six
from lxml import etree
import re
import os.path

from mpas_analysis.shared.containers import ReadOnlyDict
from mpas_analysis.shared.io.utility import paths
from mpas_analysis.shared.timekeeping.utility import string_to_datetime


def convert_namelist_to_dict(fname, readonly=True):
    """
    Converts a namelist file to key-value pairs in dictionary.

    Parameters
    ----------
    fname : str
        The file name of the namelist

    readonly : bool, optional
        Should the resulting dictionary read-only?
    """
    # Authors
    # -------
    # Phillip J Wolfram

    # form dictionary
    nml = dict()

    regex = re.compile(r"^\s*(.*?)\s*=\s*['\"]*(.*?)['\"]*\s*\n")
    with open(fname) as f:
        for line in f:
            match = regex.findall(line)
            if len(match) > 0:
                # assumes that there is only one match per line
                nml[match[0][0].lower()] = match[0][1]
    if readonly:
        nml = ReadOnlyDict(nml)

    return nml


class NameList:
    """
    Class for fortran manipulation of namelist files, provides
    read and write functionality
    """
    # Authors
    # -------
    # Phillip Wolfram, Xylar Asay-Davis

    # constructor
    def __init__(self, fname, path=None):
        """
        Parse the namelist file

        Parameters
        ----------
        fname : str
            The file name of the namelist file

        path : str, optional
            If ``fname`` contains a relative path, ``fname`` is
            relative to ``path``, rather than the current working directory
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

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
        """
        Accessor for dot noation, e.g., nml.field

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : str
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        return self.nml[key]

    # provide accessor for dictionary notation (returns string)
    def __getitem__(self, key):
        """
        Accessor for bracket noation, e.g., nml['field']

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : str
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        return self.nml[key]

    # provide accessors for get, getint, getfloat, getbool with appropriate
    # casting for comparable behavior with config files #{{{
    def get(self, key):
        """
        Get the value associated with a given key

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : str
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        return self.nml[key]

    def getint(self, key):
        """
        Get the integer value associated with a given key

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : int
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        return int(self.nml[key])

    def getfloat(self, key):
        """
        Get the float value associated with a given key

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : float
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        return float(self.nml[key])

    def getbool(self, key):
        """
        Get the boolean value associated with a given key

        Parameters
        ----------
        key : str
            The key to get a value for

        Returns
        -------
        value : bool
            The value associated with ``key``
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        if 'True' in self.nml[key] or 'true' in self.nml[key]:
            return True
        else:
            return False

    def find_option(self, possibleOptions):
        """
        If one (or more) of the names in ``possibleOptions`` is an option in
        this namelist file, returns the first match.

        Parameters
        ----------
        possibleOptions: list of str
            A list of options to search for

        Returns
        -------
        optionName : str
            The name of an option from possibleOptions occurring in the
            namelist file

        Raises
        ------
        ValueError
            If no match is found.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        for optionName in possibleOptions:
            if optionName in self.nml.keys():
                return optionName

        raise ValueError('None of the possible options {} found in namelist '
                         'file {}.'.format(possibleOptions, self.fname))

    # }}}


class StreamsFile:
    """
    Class to read in streams configuration file, provdies
    read and write functionality
    """
    # Authors
    # -------
    # Phillip Wolfram, Xylar Asay-Davis

    def __init__(self, fname, streamsdir=None):
        """
        Parse the streams file.

        Parameters
        ----------
        fname : str
            The file name the stream file

        streamsdir : str, optional
            The base path to both the output streams data and the sreams file
            (the latter only if ``fname`` is a relative path).
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

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
        """
        Get the value of the given attribute in the given stream

        Parameters
        ----------
        streamname : str
            The name of the stream

        attribname : str
            The name of the attribute within the stream


        Returns
        -------
        value : str
            The value associated with the attribute, or ``None`` if the
            attribute was not found
        """
        # Authors
        # -------
        # Phillip Wolfram, Xylar Asay-Davis

        for stream in self.root:
            # assumes streamname is unique in XML
            if stream.get('name') == streamname:
                return stream.get(attribname)
        return None

    def read_datetime_template(self, streamname):
        """
        Get the value of the given attribute in the given stream

        Parameters
        ----------
        streamname : str
            The name of the stream

        Returns
        -------
        value : str
            The template for file names from this stream in a format accepted
            by ``datetime.strptime``.  This is useful for parsing the date
            from a given file name.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        template = self.read(streamname, 'filename_template')
        replacements = {'$Y': '%Y',
                        '$M': '%m',
                        '$D': '%d',
                        '$S': '00000',  # datetime doesn't handle seconds alone
                        '$h': '%H',
                        '$m': '%M',
                        '$s': '%S'}

        for old in replacements:
            template = template.replace(old, replacements[old])

        return template

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

        calendar : {'gregorian', 'gregorian_noleap'}, optional
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
        """
        # Authors
        # -------
        # Xylar Asay-Davis

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

        if startDate is not None:
            # read one extra file before the start date to be on the safe side
            if isinstance(startDate, six.string_types):
                startDate = string_to_datetime(startDate)

        if endDate is not None:
            # read one extra file after the end date to be on the safe side
            if isinstance(endDate, six.string_types):
                endDate = string_to_datetime(endDate)

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
        Does the stream file have the given stream?

        Returns True if the streams file has a stream with the given
        streamName, otherwise returns False.

        Parameters
        ----------
        streamName : str
            The name of the stream

        Returns
        -------
        streamFound : bool
            ``True`` if the stream was found in the stream file, ``False``
            otherwise
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        for stream in self.root:
            # assumes streamname is unique in XML
            if stream.get('name') == streamName:
                return True
        return False

    def find_stream(self, possibleStreams):
        """
        If one (or more) of the names in ``possibleStreams`` is an stream in
        this streams file, returns the first match.

        Parameters
        ----------
        possibleStreams : list of str
            A list of streams to search for

        Returns
        -------
        streamName : str
            The name of an stream from possibleOptions occurring in the
            streams file

        Raises
        ------
        ValueError
            If no match is found.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        for streamName in possibleStreams:
            if self.has_stream(streamName):
                return streamName

        raise ValueError('None of the possible streams {} found in streams '
                         'file {}.'.format(possibleStreams, self.fname))

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
