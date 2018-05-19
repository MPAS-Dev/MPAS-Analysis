#!/usr/bin/env python
"""
This script is used to convert an xml-based table into an rst table for
python documentation and pythonic parsing.

Phillip J. Wolfram
Xylar Asay-Davis
04/16/2018
"""

import xml.etree.ElementTree as ET
import tabulate
import argparse
import re
import os



def markdown_links(data, footer):
    urlscmd = re.findall(r"\[.*?\]\n*\(.*?\)", data)
    urls = re.findall(r"\[.*?\]\n*\((.*?)\)", data)
    linknames = re.findall(r"\[(.*?)\]\n*\(.*?\)", data)

    for alinkname, aurl, aurlscmd in zip(linknames, urls, urlscmd):
        data = data.replace(aurlscmd, '`' + alinkname + '`_')
        footer += ".. _`" + alinkname + "`: " + aurl + "\n"

    return data, footer


def spurious_newline_whitespace(data, _):
    data = '\n'.join([string.strip() for string in data.split('\n')])
    return data, _


def cleanup(linedata, footer):
    cleanups = [spurious_newline_whitespace, markdown_links]
    for acleanup in cleanups:
        linedata, footer = acleanup(linedata, footer)
    return linedata, footer


def add_task_links(data):
    lines = data.split('\n')
    outLines = []
    for line in lines:
        line = line.strip()
        if line[0:2] == '- ':
            line = '- :ref:`task_{}`'.format(line[2:])
        else:
            line = ':ref:`task_{}`'.format(line)
        outLines.append(line)

    data = '\n'.join(outLines)
    return data


def build_rst_table_from_xml(xmlfile, rstfile, component):

    # open xml file for reading
    xml = ET.parse(xmlfile)

    # open rst file for writing
    rst = open(rstfile, 'w')

    # get headers and build list to write entries
    headers = xml.getroot().attrib['headers'].replace(' ', '').split(',')
    headernames = [aname.strip() for aname in
                   xml.getroot().attrib['headernames'].split(',')]
    data = []
    footer = ''
    for entry in xml.findall('observation'):
        if (component != 'all' and
                entry.findall('component')[0].text.strip() != component):
            continue
        line = []
        for aheader in headers:
            linedata = entry.findall(aheader)[0].text.strip()
            if aheader == 'name':
                # add a link to the associated page, if available
                nameInDocs = entry.findall('nameInDocs')
                if len(nameInDocs) > 0:
                    nameInDocs = nameInDocs[0].text.strip()
                    linedata = ':ref:`{}`'.format(nameInDocs)
            if aheader == 'tasks':
                linedata = add_task_links(linedata)
            linedata, footer = cleanup(linedata, footer)
            line.append(linedata)
        data.append(line)

    rst.writelines(tabulate.tabulate(data, headernames, tablefmt='rst') + '\n')
    rst.write('\n')
    rst.write(footer)
    rst.write('\n')

    rst.close()


def build_obs_pages_from_xml(xmlfile):

    # open xml file for reading
    xml = ET.parse(xmlfile)

    titles = {'description': 'Description',
              'source': 'Source',
              'releasePolicy': 'Release Policy',
              'references': 'References',
              'tasks': 'MPAS-Analysis Tasks'}

    path = 'obs'
    try:
        os.makedirs(path)
    except OSError:
        pass

    all_obs = open('all_obs.rst', 'w')

    all_obs.write('.. toctree::\n')
    all_obs.write('   :maxdepth: 1\n\n')

    urlBase = 'https://web.lcrc.anl.gov/public/e3sm/diagnostics/observations'

    for entry in xml.findall('observation'):

        nameInDocs = entry.findall('nameInDocs')
        if len(nameInDocs) == 0:
            # we're not making a page for this one
            continue

        footer = ''
        nameInDocs = nameInDocs[0].text.strip()

        all_obs.write('   {}/{}.rst\n'.format(path, nameInDocs))

        name = entry.findall('name')[0].text.strip()

        rst = open('{}/{}.rst'.format(path, nameInDocs), 'w')

        rst.write('.. _{}:\n\n'.format(nameInDocs))
        rst.write('{}\n'.format(name))
        rst.write('='*len(name))
        rst.write('\n\n')

        for tag in titles:
            title = titles[tag]
            rst.write('{}\n'.format(title))
            rst.write('{}\n'.format('-'*len(title)))
            text = entry.findall(tag)[0].text.strip()
            if tag == 'tasks':
                text = add_task_links(text)
            if tag == 'references':
                # add a link to the bibtex file
                subdirectory = entry.findall('subdirectory')[0].text.strip()
                text = '{}\n\n[bibtex file]({}/{}/obs.bib)'.format(
                        text, urlBase, subdirectory)
            text, footer = cleanup(text, footer)
            rst.write('{}\n\n'.format(text))
        rst.write(footer)
        rst.write('\n')
        rst.close()
    all_obs.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
            description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-x", "--xml_table", dest="xml_table",
                        help="Path to file containing "
                             "xml description of table",
                        metavar="FILE", required=True)
    parser.add_argument("-r", "--rst_table", dest="rst_table",
                        help="Path to file containing rst description "
                             "of table for output",
                        metavar="FILE", required=True)
    parser.add_argument("-c", "--component", dest="component",
                        help="Component for parsing of table: "
                             "'landice', 'ocean', 'seaice', or 'all'",
                        metavar="STRING", default="all")

    args = parser.parse_args()

    build_rst_table_from_xml(args.xml_table, args.rst_table, args.component)
