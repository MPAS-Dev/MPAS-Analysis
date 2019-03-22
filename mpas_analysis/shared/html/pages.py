# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2019 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2019 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2019 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/master/LICENSE

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import pkg_resources
from os import makedirs
from shutil import copyfile
from lxml import etree
from collections import OrderedDict
import subprocess
import os

import mpas_analysis
from mpas_analysis.shared.io.utility import build_config_full_path


def generate_html(config, analyses, controlConfig=None):  # {{{
    """
    Generates webpages for diplaying the plots from each analysis task

    Parameters
    ----------
    config : ``MpasAnalysisConfigParser``
        Config options

    analysis : ``OrderedDict`` of ``AnalysisTask`` objects
        the analysis tasks that generated the plots to include in the webpages.
        The ``list_xml_files()`` method will be called on each task to get
        the list of files to include on the webpage for the associated
        component.

    controlConfig : ``MpasAnalysisConfigParser``, optional
        Config options for a control run

    """
    # Authors
    # -------
    # Xylar Asay-Davis

    generateHTML = config.getboolean('html', 'generate')
    if not generateHTML:
        return

    print("Generating webpage for viewing results...")

    page = MainPage(config, controlConfig)

    components = OrderedDict()

    # add images from each analysis task, creating ga dictionary of components
    missingCount = 0
    for analysisTask in analyses.values():
        for fileName in analysisTask.xmlFileNames:
            try:
                ComponentPage.add_image(fileName, config, components,
                                        controlConfig)
            except IOError:
                print('  missing file {}'.format(fileName))
                missingCount += 1

    if missingCount > 0:
        print('Warning: {} XML files were missing and the analysis website'
              ' will be incomplete.'.format(missingCount))
    # generate the page for each component and add the component to the main
    # page
    for componentName, component in components.items():
        component.generate()

        firstImageFileName = component.get_first_image()

        page.add_component(componentName, component.subdirectory,
                           firstImageFileName)

    page.generate()

    print("Done.")

    # }}}


class MainPage(object):
    """
    Describes a main webpage containg one or more pages for components

    Attributes
    ----------
    config : ``MpasAnalysisConfigParser``
        Config options

    controlConfig : ``MpasAnalysisConfigParser``
        Config options for a control run

    pageTemplate, componentTemplate : str
        The contents of templates used to construct the page

    components : OrederdDict of dict
        Each component has a name, subdirectory and image name used to find
        the appropriate thumbnail.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, controlConfig=None):
        """
        Create a MainPage object, reading in the templates

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Config options

        controlConfig : ``MpasAnalysisConfigParser``, optional
            Config options for a control run
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.config = config
        self.controlConfig = controlConfig

        # get template text
        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/main_page.html")

        with open(fileName, 'r') as templateFile:
            self.pageTemplate = templateFile.read()

        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/main_component.html")
        with open(fileName, 'r') as templateFile:
            self.componentTemplate = templateFile.read()

        # start with no components
        self.components = OrderedDict()

    def add_component(self, name, subdirectory, imageFileName):
        """
        Create a MainPage object, reading in the templates

        Parameters
        ----------
        name : str
            The name of the component as it should appear in the list of
            components, at the top of the component webpage and in the page
            title (e.g "Sea Ice" as opposed to "sea_ice" or "seaIce")

        subdirecory : str
            The subdirectory for the component's webpage

        imageFileName : str
            The name of an image file (without path) that will be used as the
            thumbnail for the gallery.  Typically, this is the first image
            from the first gallery.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.components[name] = {'subdirectory': subdirectory,
                                 'imageFileName': imageFileName}

    def generate(self):
        """
        Generate the webpage from templates and components, and write it out to
        the HTML directory.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        runName = self.config.get('runs', 'mainRunName')

        if self.controlConfig is None:
            controlRunText = ''
        else:
            controlRunText = '<br> Control: {}'.format(
                self.controlConfig.get('runs', 'mainRunName'))

        componentsText = ''

        for componentName, componentDict in self.components.items():
            subdirectory = componentDict['subdirectory']
            imageFileName = componentDict['imageFileName']
            replacements = {'@componentDir': subdirectory,
                            '@componentName': componentName,
                            '@firstImage': imageFileName}

            # substitute entries in the template and add the component to
            # the text describing all components
            componentsText = componentsText + \
                _replace_tempate_text(self.componentTemplate, replacements)

        githash = _get_git_hash()
        if githash is None:
            githash = ''
        else:
            githash = 'Git Hash: {}'.format(githash)

        replacements = {'@runName': runName,
                        '@controlRunText': controlRunText,
                        '@components': componentsText,
                        '@version': mpas_analysis.__version__,
                        '@gitHash': githash}

        pageText = _replace_tempate_text(self.pageTemplate, replacements)

        htmlBaseDirectory = build_config_full_path(self.config, 'output',
                                                   'htmlSubdirectory')

        for subdir in ['css', 'js']:
            try:
                makedirs('{}/{}'.format(htmlBaseDirectory, subdir))
            except OSError:
                pass

        outFileName = '{}/index.html'.format(htmlBaseDirectory)

        with open(outFileName, mode='w') as mainFile:
            mainFile.write(
                pageText.encode('ascii',
                                'xmlcharrefreplace').decode('ascii'))

        # copy the css and js files as well as general images
        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/style.css")
        copyfile(fileName, '{}/css/style.css'.format(htmlBaseDirectory))

        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/index.js")
        copyfile(fileName, '{}/js/index.js'.format(htmlBaseDirectory))

        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/mpas_logo.png")
        copyfile(fileName, '{}/mpas_logo.png'.format(htmlBaseDirectory))

        fileName = \
            pkg_resources.resource_filename(__name__,
                                            "templates/config.png")
        copyfile(fileName, '{}/config.png'.format(htmlBaseDirectory))

        with open('{}/config.{}'.format(htmlBaseDirectory, runName), 'w') \
                as configFile:
            self.config.write(configFile)


class ComponentPage(object):
    """
    Describes a component with one or more gallery groups, each with one or
    more galleries, and a list of "quick links" to the gallery groups

    Attributes
    ----------
    config : ``MpasAnalysisConfigParser``
        Config options

    controlConfig : ``MpasAnalysisConfigParser``
        Config options for a control run

    name : str
        The name of the component as it should appear in the list of
        components, at the top of the component webpage and in the page
        title (e.g "Sea Ice" as opposed to "sea_ice" or "seaIce")

    subdirecory : str
        The subdirectory for the component's webpage

    templates : OrderedDict of str
        The contents of templates used to construct the page

    groups : tree of OrederdDict
        A tree of information describing the the gallery groups in the page,
        the galleries in each group and the images in each gallery.
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    def __init__(self, config, name, subdirectory, controlConfig=None):
        """
        Create a ComponentPage object, reading in the templates

        Parameters
        ----------
        config : ``MpasAnalysisConfigParser``
            Config options

        name : str
            The name of the component as it should appear in the list of
            components, at the top of the component webpage and in the page
            title (e.g "Sea Ice" as opposed to "sea_ice" or "seaIce")

        subdirecory : str
            The subdirectory for the component's webpage

        controlConfig : ``MpasAnalysisConfigParser``, optional
            Config options for a control run
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        self.config = config
        self.controlConfig = controlConfig
        self.name = name
        self.subdirectory = subdirectory

        htmlBaseDirectory = build_config_full_path(self.config, 'output',
                                                   'htmlSubdirectory')
        self.directory = '{}/{}'.format(htmlBaseDirectory, self.subdirectory)

        self.templates = OrderedDict()

        for templateName in ['page', 'quicklink', 'group', 'gallery', 'image',
                             'subtitle']:

            # get template text
            fileName = pkg_resources.resource_filename(
                __name__,
                "templates/component_{}.html".format(templateName))

            with open(fileName, 'r') as templateFile:
                self.templates[templateName] = templateFile.read()

        # start with no groups
        self.groups = OrderedDict()

    @staticmethod
    def add_image(xmlFileName, config, components, controlConfig=None):
        """
        Add the image to the appropriate component.  Note: this is a static
        method because we do not know which component to add the image to
        until we have read the XML file.

        Parameters
        ----------
        xmlFileName : str
            The full path to the XML file describing the image to be added

        config : ``MpasAnalysisConfigParser`` object
            contains config options

        components : OrederdDict of dict
            A dictionary of components to which the image will be added.  If
            the appropriate component is not yet in the dictionary, it will
            be added. ``components`` should be viewed as an input and output
            parameter, since it is modified by this function.

        controlConfig : ``MpasAnalysisConfigParser``, optional
            Config options for a control run
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        xmlRoot = etree.parse(xmlFileName).getroot()

        componentName = ComponentPage._get_required_xml_text(xmlRoot,
                                                             'componentName',
                                                             xmlFileName)
        componentSubdirectory = ComponentPage._get_required_xml_text(
            xmlRoot, 'componentSubdirectory', xmlFileName)

        imageFileName = ComponentPage._get_required_xml_text(xmlRoot,
                                                             'imageFileName',
                                                             xmlFileName)
        groupName = ComponentPage._get_required_xml_text(xmlRoot,
                                                         'galleryGroup',
                                                         xmlFileName)
        groupLink = ComponentPage._get_required_xml_text(xmlRoot,
                                                         'groupLink',
                                                         xmlFileName)

        if componentName not in components:
            components[componentName] = ComponentPage(config, componentName,
                                                      componentSubdirectory,
                                                      controlConfig)

        component = components[componentName]

        if groupName not in component.groups:
            component.groups[groupName] = {'galleries': OrderedDict(),
                                           'link': groupLink}
            group = component.groups[groupName]
            node = xmlRoot.find('groupSubtitle')
            if node is not None:
                group['subtitle'] = node.text

        node = xmlRoot.find('gallery')
        if node is None:
            galleryName = 'None'
        else:
            galleryName = node.text

        galleries = component.groups[groupName]['galleries']
        if galleryName not in galleries:
            galleries[galleryName] = {'images': OrderedDict()}

        images = galleries[galleryName]['images']
        if imageFileName in images:
            raise ValueError('image {} already added to component page '
                             '{}'.format(imageFileName, componentName))

        images[imageFileName] = OrderedDict()
        image = images[imageFileName]
        for tag in ['thumbnailDescription', 'imageDescription',
                    'imageCaption', 'imageSize', 'orientation']:
            node = xmlRoot.find(tag)
            if node is None or node.text is None:
                image[tag] = ''
            else:
                image[tag] = node.text

    def generate(self):
        """
        Generate the webpage from templates and groups, and write it out to
        the HTML directory.
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        runName = self.config.get('runs', 'mainRunName')

        if self.controlConfig is None:
            controlRunText = ''
        else:
            controlRunText = '<br> Control: {}'.format(
                self.controlConfig.get('runs', 'mainRunName'))

        quickLinkText = ''
        galleriesText = ''
        for groupName, groupDict in self.groups.items():
            quickLinkText = quickLinkText + \
                self._generate_quick_link_text(groupName, groupDict)

            galleriesText = galleriesText + \
                self._generate_group_text(groupName, groupDict)

        replacements = {'@runName': runName,
                        '@controlRunText': controlRunText,
                        '@componentName': self.name,
                        '@quickLinks': quickLinkText,
                        '@galleries': galleriesText}

        pageText = _replace_tempate_text(self.templates['page'], replacements)

        outFileName = '{}/index.html'.format(self.directory)

        with open(outFileName, mode='w') as componentFile:
            componentFile.write(
                pageText.encode('ascii',
                                'xmlcharrefreplace').decode('ascii'))

    def get_first_image(self):
        """
        Find the first image in the first gallery in this component (typically
        to use it as a thumbnail)

        Returns
        -------
        firstImageFilename : str
            The name (with out path) of the first image in the first gallery
        """
        # Authors
        # -------
        # Xylar Asay-Davis

        # get the first image name
        firstGroup = next(iter(self.groups.values()))
        firstGallery = next(iter(firstGroup['galleries'].values()))
        firstImageFileName = next(iter(firstGallery['images']))
        return firstImageFileName

    @staticmethod
    def _get_required_xml_text(root, tag, fileName):
        """read the value associated with a required tag from the XML root"""
        node = root.find(tag)
        if node is None:
            raise IOError('image descriptor file {} is missing a required'
                          '{} entry'.format(fileName, tag))
        return node.text

    def _generate_image_text(self, imageFileName, imageDict):
        """fill in the template for a given image with the desired content"""
        replacements = {'@imageFileName': imageFileName}
        for tag in ['imageSize', 'imageDescription', 'imageCaption',
                    'thumbnailDescription', 'orientation']:
            replacements['@{}'.format(tag)] = imageDict[tag]

        imageText = _replace_tempate_text(self.templates['image'],
                                          replacements)
        return imageText

    def _generate_gallery_text(self, galleryName, images):
        """fill in the template for a given gallery with the desired content"""
        imagesText = ''
        for imageFileName, imageDict in images.items():
            imagesText = imagesText + \
                self._generate_image_text(imageFileName, imageDict)

        if galleryName == 'None':
            galleryTitle = ''
        else:
            galleryTitle = self._generate_subtitle_text(galleryName)

        replacements = {'@galleryTitle': galleryTitle,
                        '@galleryImages': imagesText}
        galleryText = _replace_tempate_text(self.templates['gallery'],
                                            replacements)
        return galleryText

    def _generate_subtitle_text(self, subtitle):
        """
        fill in the template for a given gallery subtitle with the desired
        content
        """
        replacements = {'@subtitle': subtitle}

        subtitleText = _replace_tempate_text(self.templates['subtitle'],
                                             replacements)
        return subtitleText

    def _generate_group_text(self, groupName, groupDict):
        """
        fill in the template for a given gallery group with the desired
        content
        """
        galleriesText = ''
        for galleryName, galleryDict in groupDict['galleries'].items():
            galleriesText = galleriesText + \
                self._generate_gallery_text(galleryName, galleryDict['images'])

        replacements = {'@analysisGroupName': groupName,
                        '@analysisGroupLink': groupDict['link'],
                        '@groupGalleries': galleriesText}

        if 'subtitle' in groupDict:
            subtitleText = self._generate_subtitle_text(groupDict['subtitle'])
        else:
            subtitleText = ''

        replacements['@groupSubtitle'] = subtitleText

        groupText = _replace_tempate_text(self.templates['group'],
                                          replacements)
        return groupText

    def _generate_quick_link_text(self, groupName, groupDict):
        """
        fill in the template for a given quick link with the desired
        content
        """

        firstGallery = next(iter(groupDict['galleries'].values()))
        firstImageFileName = next(iter(firstGallery['images']))

        replacements = {'@analysisGroupName': groupName,
                        '@analysisGroupLink': groupDict['link'],
                        '@imageFileName': firstImageFileName}

        quickLinkText = _replace_tempate_text(self.templates['quicklink'],
                                              replacements)

        return quickLinkText


def _replace_tempate_text(template, replacements):
    """
    replace substrings in a given template based on a dictionary of
    replacements.
    """
    output = template
    for src, target in replacements.items():
        output = output.replace(src, target)
    return output


def _get_git_hash():
    """
    get the hashtag of the git commit (if any) that MPAS-Analysis is being run
    from
    """
    with open(os.devnull, 'w') as devnull:
        try:
            githash = subprocess.check_output(['git', 'log',
                                               '--pretty=format:"%h"',
                                               '-n', '1'],
                                              stderr=devnull)
        except subprocess.CalledProcessError:
            return None

    githash = githash.decode('utf-8').strip('\n').replace('"', '')
    return githash
