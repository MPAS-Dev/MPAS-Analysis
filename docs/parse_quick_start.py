#!/usr/bin/env python
"""
A script for converting the README.md to a quick-start guide for inclusion
in the documentation

Xylar Asay-Davis
11/06/2018
"""

from m2r import convert


def build_quick_start():

    replace = {'# MPAS-Analysis': '# Quick Start Guide\n',
               '[![Build Status]': '',
               '[![Documentation Status]': '',
               '## Documentation': '',
               '[http://mpas-analysis.readthedocs.io]': '',
               '![sea surface temperature](docs/_static/sst_example.png)':
                   '![sea surface temperature](_static/sst_example.png)\n'}
    outContent = ''
    with open('../README.md', 'r') as inFile:
        for line in inFile.readlines():
            for replaceString in replace:
                if replaceString in line:
                    line = replace[replaceString]
                    break
            outContent = outContent + line

    outContent = convert(outContent)

    with open('quick_start.rst', 'w') as outFile:
        outFile.write('.. _quick_start:\n\n')
        outFile.write(outContent)
