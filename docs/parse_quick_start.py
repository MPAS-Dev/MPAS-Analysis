#!/usr/bin/env python
"""
A script for converting the README.md to a quick-start guide for inclusion
in the documentation

Xylar Asay-Davis
05/15/2018
"""


def build_quick_start():

    replace = {'# MPAS-Analysis': '# Quick Start Guide\n',
               '[![Build Status]': '',
               '[![Documentation Status]': '',
               '## Documentation': '',
               '[http://mpas-analysis.readthedocs.io]': '',
               '![sea surface temperature](docs/_static/sst_example.png)':
                   '![sea surface temperature](_static/sst_example.png)\n'}
    outContent = []
    with open('../README.md', 'r') as inFile:
        with open('quick_start.md', 'w') as outFile:
            outFile.write(
                '```eval_rst\n'
                '.. _quick_start:\n'
                '```\n')
            for line in inFile.readlines():
                for replaceString in replace:
                    if replaceString in line:
                        line = replace[replaceString]
                        break
                outContent.append(line)
            outFile.writelines(outContent)
