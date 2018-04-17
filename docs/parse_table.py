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
import re

def markdown_links(data, footer):
    urlscmd = re.findall(r"\[.*?\]\n*\(.*?\)", data)
    urls = re.findall(r"\[.*?\]\n*\((.*?)\)", data)
    linknames = re.findall(r"\[(.*?)\]\n*\(.*?\)", data)

    for alinkname, aurl, aurlscmd in zip(linknames, urls, urlscmd):
        data = data.replace(aurlscmd, '`' + alinkname +'`_')
        footer += ".. _`" + alinkname + "`: " + aurl + "\n"

    return data, footer

def spurious_newline_whitespace(data, _):
    whitespace = re.findall('\n\s*', data)
    if len(whitespace) > 0:
        astr = min(whitespace)
        data = data.replace(astr, "\n")
    return data, _

def cleanup(linedata, footer):
    cleanups = [spurious_newline_whitespace, markdown_links]
    for acleanup in cleanups:
        linedata, footer = acleanup(linedata, footer)
    return linedata, footer

def build_rst_table_from_xml(xmlfile, rstfile, component):

    # open xml file for reading
    xml = ET.parse(xmlfile)

    # open rst file for writing
    rst = open(rstfile, 'w')

    # get headers and build list to write entries
    headers = xml.getroot().attrib['headers'].replace(' ','').split(',')
    headernames = [aname.strip() for aname in xml.getroot().attrib['headernames'].split(',')]
    data = []
    footer = '\n'
    for entry in xml.findall('observation'):
        if component != 'all' and entry.findall('component')[0].text.strip() != component:
            continue
        line = []
        for aheader in headers:
            linedata = entry.findall(aheader)[0].text.strip()
            linedata, footer = cleanup(linedata, footer)
            line.append(linedata)
        data.append(line)

    rst.write(tabulate.tabulate(data, headernames, tablefmt='rst'))
    rst.write(footer)

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
    parser.add_argument("-c", "--component", dest="component",
                        help="Component for parsing of table, 'landice', 'ocean', "
                             "'seaice', or 'all'", metavar="STRING", default="all")

    args = parser.parse_args()

    build_rst_table_from_xml(args.xml_table, args.rst_table, args.component)
