#!/usr/bin/env python
"""
This script is used to convert an xml-based table into an rst table for
python documentation and pythonic parsing.

Phillip J. Wolfram
04/16/2018
"""

import xml.etree.ElementTree as ET
import tabulate
import argparse

def build_rst_table_from_xml(xmlfile, rstfile):

    # open xml file for reading
    xml = ET.parse(xmlfile)

    # open rst file for writing
    rst = open(rstfile, 'w')

    # get headers and build list to write entries
    headers = xml.getroot().attrib['headers'].replace(' ','').split(',')
    headernames = [aname.strip() for aname in xml.getroot().attrib['headernames'].split(',')]
    data = []
    for entry in xml.findall('aobs'):
        line = []
        for aheader in headers:
            line.append(entry.findall(aheader)[0].text.strip().replace('\n', ' '))
        data.append(line)

    rst.write(tabulate.tabulate(data, headernames, tablefmt='rst'))

    rst.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-x", "--xml_table", dest="xml_table",
                        help="Path to file containing xml description of table",
                        metavar="FILE", required=True)
    parser.add_argument("-r", "--rst_table", dest="rst_table",
                        help="Path to file containing rst description of table for output",
                        metavar="FILE", required=True)

    args = parser.parse_args()

    build_rst_table_from_xml(args.xml_table, args.rst_table)
