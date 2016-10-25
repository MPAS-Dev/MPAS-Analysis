#!/usr/bin/env python
"""
    Module of classes / routines to manipulate fortran namelist and streams
    files.

    Phillip Wolfram
    LANL
    04/24/2014
"""

from subprocess import check_output, call
from lxml import etree
import re
from ..containers import ReadOnlyDict

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
    """

    # constructor
    def __init__(self, fname):
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
    """

    def __init__(self, fname):
        self.fname = fname
        self.xmlfile = etree.parse(fname)
        self.root = self.xmlfile.getroot()

    def read(self, streamname, attribname):
        """ name is a list of name entries terminanting in some value
        """
        for stream in self.root:
            # assumes streamname is unique in XML
            if stream.get('name') == streamname:
                return stream.get(attribname)


# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
