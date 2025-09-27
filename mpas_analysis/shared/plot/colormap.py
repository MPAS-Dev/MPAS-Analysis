# This software is open source software available under the BSD-3 license.
#
# Copyright (c) 2022 Triad National Security, LLC. All rights reserved.
# Copyright (c) 2022 Lawrence Livermore National Security, LLC. All rights
# reserved.
# Copyright (c) 2022 UT-Battelle, LLC. All rights reserved.
#
# Additional copyright and license information can be found in the LICENSE file
# distributed with this code, or at
# https://raw.githubusercontent.com/MPAS-Dev/MPAS-Analysis/main/LICENSE
"""
Utilities for handling color maps and color bars
"""
# Authors
# -------
# Xylar Asay-Davis, Milena Veneziani, Luke Van Roekel, Greg Streletz

import matplotlib.pyplot as plt
import matplotlib.colors as cols
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib
from xml.etree import ElementTree
import configparser
import cmocean
import importlib.resources


def setup_colormap(config, configSectionName, suffix=''):
    """
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    Returns
    -------
    colormapDict : dict
        A dictionary of colormap information.

        'colormap' specifies the name of the new colormap

        'norm' is a matplotlib norm object used to normalize the colormap

        'levels' is an array of contour levels or ``None`` if not using indexed
        color map

        'ticks' is an array of values where ticks should be placed

        'contours' is an array of contour values to plot or ``None`` if none
        have been specified

        'lineWidth' is the width of contour lines or ``None`` if not specified

        'lineColor' is the color of contour lines or ``None`` if not specified

        `arrowWidth` is the width of the arrows or ``None`` if not specified

        `arrowSpacing` is the spacing between arrows or ``None`` if not
        specified
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    register_custom_colormaps()

    option = f'colormapType{suffix}'
    if config.has_option(configSectionName, option):
        colormapType = config.get(configSectionName, option)
        if colormapType == 'indexed':
            (colormap, norm, levels, ticks) = _setup_indexed_colormap(
                config, configSectionName, suffix=suffix)
        elif colormapType == 'continuous':
            (colormap, norm, ticks) = _setup_colormap_and_norm(
                config, configSectionName, suffix=suffix)
            levels = None
        else:
            raise ValueError(f'config section {configSectionName} option '
                             f'{option} is not "indexed" or "continuous"')
    else:
        colormap = None
        norm = None
        levels = None
        ticks = None

    contours = None
    lineWidth = None
    lineColor = None
    arrowWidth = None
    arrowSpacing = None

    option = f'contourLevels{suffix}'
    if config.has_option(configSectionName, option):
        contours = config.getnumpy(configSectionName, option)
        if isinstance(contours, str) and contours == 'none':
            contours = None

    option = f'contourThickness{suffix}'
    if config.has_option(configSectionName, option):
        lineWidth = config.getfloat(configSectionName, option)

    option = f'contourColor{suffix}'
    if config.has_option(configSectionName, option):
        lineColor = config.get(configSectionName, option)

    option = f'arrowWidth{suffix}'
    if config.has_option(configSectionName, option):
        arrowWidth = config.getexpression(configSectionName, option)
        if isinstance(arrowWidth, str) and arrowWidth == 'None':
            arrowWidth = None
    option = f'arrowSpacing{suffix}'
    if config.has_option(configSectionName, option):
        arrowSpacing = config.getexpression(configSectionName, option)
        if isinstance(arrowSpacing, str) and arrowSpacing == 'None':
            arrowSpacing = None

    return {'colormap': colormap, 'norm': norm, 'levels': levels,
            'ticks': ticks, 'contours': contours, 'lineWidth': lineWidth,
            'lineColor': lineColor, 'arrowWidth': arrowWidth,
            'arrowSpacing': arrowSpacing}


def register_custom_colormaps():
    name = 'ferret'
    backgroundColor = (0.9, 0.9, 0.9)

    red = np.array([[0, 0.6],
                    [0.15, 1],
                    [0.35, 1],
                    [0.65, 0],
                    [0.8, 0],
                    [1, 0.75]])

    green = np.array([[0, 0],
                      [0.1, 0],
                      [0.35, 1],
                      [1, 0]])

    blue = np.array([[0, 0],
                     [0.5, 0],
                     [0.9, 0.9],
                     [1, 0.9]])

    colorCount = 21
    colorList = np.ones((colorCount, 4), float)
    colorList[:, 0] = np.interp(np.linspace(0, 1, colorCount),
                                red[:, 0], red[:, 1])
    colorList[:, 1] = np.interp(np.linspace(0, 1, colorCount),
                                green[:, 0], green[:, 1])
    colorList[:, 2] = np.interp(np.linspace(0, 1, colorCount),
                                blue[:, 0], blue[:, 1])
    colorList = colorList[::-1, :]

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    colorMap.set_bad(backgroundColor)
    _register_colormap_and_reverse(name, colorMap)

    name = 'erdc_iceFire_H'

    colorArray = np.array([
        [-1, 4.05432e-07, 0, 5.90122e-06],
        [-0.87451, 0, 0.120401, 0.302675],
        [-0.74902, 0, 0.216583, 0.524574],
        [-0.623529, 0.0552475, 0.345025, 0.6595],
        [-0.498039, 0.128047, 0.492588, 0.720288],
        [-0.372549, 0.188955, 0.641309, 0.792092],
        [-0.247059, 0.327673, 0.784935, 0.873434],
        [-0.121569, 0.60824, 0.892164, 0.935547],
        [0.00392157, 0.881371, 0.912178, 0.818099],
        [0.129412, 0.951407, 0.835621, 0.449279],
        [0.254902, 0.904481, 0.690489, 0],
        [0.380392, 0.85407, 0.510864, 0],
        [0.505882, 0.777093, 0.33018, 0.00088199],
        [0.631373, 0.672862, 0.139087, 0.00269398],
        [0.756863, 0.508815, 0, 0],
        [0.882353, 0.299417, 0.000366289, 0.000547829],
        [1, 0.0157519, 0.00332021, 4.55569e-08]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex + 1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    _register_colormap_and_reverse(name, colorMap)

    name = 'erdc_iceFire_L'

    colorArray = np.array([
        [-1, 0.870485, 0.913768, 0.832905],
        [-0.87451, 0.586919, 0.887865, 0.934003],
        [-0.74902, 0.31583, 0.776442, 0.867858],
        [-0.623529, 0.18302, 0.632034, 0.787722],
        [-0.498039, 0.117909, 0.484134, 0.713825],
        [-0.372549, 0.0507239, 0.335979, 0.654741],
        [-0.247059, 0, 0.209874, 0.511832],
        [-0.121569, 0, 0.114689, 0.28935],
        [0.00392157, 0.0157519, 0.00332021, 4.55569e-08],
        [0.129412, 0.312914, 0, 0],
        [0.254902, 0.520865, 0, 0],
        [0.380392, 0.680105, 0.15255, 0.0025996],
        [0.505882, 0.785109, 0.339479, 0.000797922],
        [0.631373, 0.857354, 0.522494, 0],
        [0.756863, 0.910974, 0.699774, 0],
        [0.882353, 0.951921, 0.842817, 0.478545],
        [1, 0.881371, 0.912178, 0.818099]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex + 1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    _register_colormap_and_reverse(name, colorMap)

    name = 'BuOr'
    colors1 = plt.cm.PuOr(np.linspace(0., 1, 256))
    colors2 = plt.cm.RdBu(np.linspace(0, 1, 256))

    # combine them and build a new colormap, just the orange from the first
    # and the blue from the second
    colorList = np.vstack((colors1[0:128, :], colors2[128:256, :]))
    # reverse the order
    colorList = colorList[::-1, :]
    colorMap = cols.LinearSegmentedColormap.from_list(name, colorList)

    _register_colormap_and_reverse(name, colorMap)

    name = 'Maximenko'
    colorArray = np.array([
        [-1, 0., 0.45882352941, 0.76470588235],
        [-0.666667, 0., 0.70196078431, 0.90588235294],
        [-0.333333, 0.3294117647, 0.87058823529, 1.],
        [0., 0.76470588235, 0.94509803921, 0.98039215686],
        [0.333333, 1., 1., 0.],
        [0.666667, 1., 0.29411764705, 0.],
        [1, 1., 0., 0.]], float)

    colorCount = 255
    colorList = np.ones((colorCount, 4), float)
    x = colorArray[:, 0]
    for cIndex in range(3):
        colorList[:, cIndex] = np.interp(
            np.linspace(-1., 1., colorCount),
            x, colorArray[:, cIndex + 1])

    colorMap = cols.LinearSegmentedColormap.from_list(
        name, colorList, N=255)

    _register_colormap_and_reverse(name, colorMap)

    # add the cmocean color maps
    map_names = list(cmocean.cm.cmapnames)
    # don't bother with gray (already exists, I think)
    map_names.pop(map_names.index('gray'))
    for map_name in map_names:
        _register_colormap_and_reverse(map_name, getattr(cmocean.cm, map_name))

    # add ScientificColourMaps7 from
    # http://www.fabiocrameri.ch/colourmaps.php
    # https://doi.org/10.5281/zenodo.5501399
    for map_name in ['acton', 'bam', 'bamako', 'bamO', 'batlow', 'batlowK',
                     'batlowW', 'berlin', 'bilbao', 'broc', 'brocO', 'buda',
                     'bukavu', 'cork', 'corkO', 'davos', 'devon', 'fes',
                     'grayC', 'hawaii', 'imola', 'lajolla', 'lapaz', 'lisbon',
                     'nuuk', 'oleron', 'oslo', 'roma', 'romaO', 'tofino',
                     'tokyo', 'turku', 'vanimo', 'vik', 'vikO']:
        package = f'mpas_analysis.shared.plot.ScientificColourMaps7.{map_name}'
        filename = f'{map_name}_PARAVIEW.xml'
        with importlib.resources.path(package,  filename) as xml_file:
            _read_xml_colormap(xml_file, map_name)

    # add SciVisColor colormaps from
    # https://sciviscolor.org/home/colormaps/
    for map_name in ['3wave-yellow-grey-blue', '3Wbgy5',
                     '4wave-grey-red-green-mgreen', '5wave-yellow-brown-blue',
                     'blue-1', 'blue-3', 'blue-6', 'blue-8', 'blue-orange-div',
                     'brown-2', 'brown-5', 'brown-8', 'green-1', 'green-4',
                     'green-7', 'green-8', 'orange-5', 'orange-6',
                     'orange-green-blue-gray', 'purple-7', 'purple-8', 'red-1',
                     'red-3', 'red-4', 'yellow-1', 'yellow-7']:
        package = 'mpas_analysis.shared.plot.SciVisColorColormaps'
        filename = f'{map_name}.xml'
        with importlib.resources.path(package,  filename) as xml_file:
            _read_xml_colormap(xml_file, map_name)

    # add SciVisColor colormaps created by hand using
    # https://sciviscolor.org/color-moves-app/
    for map_name in ['3wave-green-red-purple', '3wave-blue-red-brown',
                     'div-one-third-blue-two-thirds-red',
                     'div-one-third-green-two-thirds-red',]:
        package = 'mpas_analysis.shared.plot.SciVisColorCustom'
        filename = f'{map_name}.xml'
        with importlib.resources.path(package,  filename) as xml_file:
            _read_xml_colormap(xml_file, map_name)

    name = 'white_cmo_deep'
    # modify cmo.deep to start at white
    colors2 = plt.get_cmap('cmo.deep')(np.linspace(0, 1, 224))
    colorCount = 32
    colors1 = np.ones((colorCount, 4), float)
    x = np.linspace(0., 1., colorCount+1)[0:-1]
    white = [1., 1., 1., 1.]
    for cIndex in range(4):
        colors1[:, cIndex] = np.interp(x, [0., 1.],
                                       [white[cIndex], colors2[0, cIndex]])

    colors = np.vstack((colors1, colors2))

    # generating a smoothly-varying LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list(name, colors)
    _register_colormap_and_reverse(name, cmap)


def _setup_colormap_and_norm(config, configSectionName, suffix=''):
    """
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    Returns
    -------
    colormap : srt
        new colormap

    norm : matplotlib.colors.Normalize
        the norm used to normalize the colormap

    ticks : array of float
        the tick marks on the colormap
    """
    # Authors
    # -------
    # Xylar Asay-Davis

    register_custom_colormaps()

    colormap = plt.get_cmap(config.get(configSectionName,
                                       f'colormapName{suffix}'))

    normType = config.get(configSectionName, f'normType{suffix}')

    kwargs = config.getexpression(configSectionName,
                                  f'normArgs{suffix}')

    if normType == 'symLog':
        norm = cols.SymLogNorm(**kwargs)
    elif normType == 'log':
        norm = cols.LogNorm(**kwargs)
    elif normType == 'linear':
        norm = cols.Normalize(**kwargs)
    else:
        raise ValueError(f'Unsupported norm type {normType} in section '
                         f'{configSectionName}')

    try:
        ticks = config.getnumpy(configSectionName, f'colorbarTicks{suffix}')
    except configparser.NoOptionError:
        ticks = None

    return colormap, norm, ticks


def _setup_indexed_colormap(config, configSectionName, suffix=''):
    """
    Set up a colormap from the registry

    Parameters
    ----------
    config : instance of ConfigParser
        the configuration, containing a [plot] section with options that
        control plotting

    configSectionName : str
        name of config section

    suffix: str, optional
        suffix of colormap related options

    Returns
    -------
    colormap : srt
        new colormap

    norm : matplotlib.colors.Normalize
        the norm used to normalize the colormap

    ticks : array of float
        the tick marks on the colormap
    """
    # Authors
    # -------
    # Xylar Asay-Davis, Milena Veneziani, Greg Streletz

    colormap = plt.get_cmap(config.get(configSectionName,
                                       f'colormapName{suffix}'))

    indices = config.getnumpy(configSectionName, f'colormapIndices{suffix}')

    try:
        levels = config.getnumpy(configSectionName, f'colorbarLevels{suffix}')
    except configparser.NoOptionError:
        levels = None

    if levels is not None:
        # set under/over values based on the first/last indices in the colormap
        underColor = colormap(indices[0])
        overColor = colormap(indices[-1])
        if len(levels) + 1 == len(indices):
            # we have 2 extra values for the under/over so make the colormap
            # without these values
            indices = indices[1:-1]
        elif len(levels) - 1 != len(indices):
            # indices list must be either one element shorter
            # or one element longer than colorbarLevels list
            raise ValueError('length mismatch between indices and '
                             'colorbarLevels')
        colormap = cols.ListedColormap(colormap(indices),
                                       f'colormapName{suffix}')
        colormap.set_under(underColor)
        colormap.set_over(overColor)

    norm = cols.BoundaryNorm(levels, colormap.N)

    try:
        ticks = config.getnumpy(configSectionName, f'colorbarTicks{suffix}')
    except configparser.NoOptionError:
        ticks = levels

    return colormap, norm, levels, ticks


def _read_xml_colormap(xmlFile, map_name):
    """Read in an XML colormap"""

    xml = ElementTree.parse(xmlFile)

    root = xml.getroot()
    colormap = root.findall('ColorMap')
    if len(colormap) > 0:
        colormap = colormap[0]
        colorDict = {'red': [], 'green': [], 'blue': []}
        for point in colormap.findall('Point'):
            x = float(point.get('x'))
            color = [float(point.get('r')), float(point.get('g')),
                     float(point.get('b'))]
            colorDict['red'].append((x, color[0], color[0]))
            colorDict['green'].append((x, color[1], color[1]))
            colorDict['blue'].append((x, color[2], color[2]))
        cmap = LinearSegmentedColormap(map_name, colorDict, 256)

        _register_colormap_and_reverse(map_name, cmap)


def _register_colormap_and_reverse(map_name, cmap):
    if map_name not in matplotlib.colormaps:
        matplotlib.colormaps.register(cmap, name=map_name)
        matplotlib.colormaps.register(cmap.reversed(), name=f'{map_name}_r')


def _plot_color_gradients():
    """from https://matplotlib.org/tutorials/colors/colormaps.html"""

    cmap_list = [m for m in plt.colormaps() if not m.endswith("_r")]

    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))

    nrows = len(cmap_list)

    fig, axes = plt.subplots(figsize=(7.2, 0.25 * nrows), nrows=nrows)
    fig.subplots_adjust(top=0.99, bottom=0.01, left=0.35, right=0.99)

    for ax, name in zip(axes, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(name))
        pos = list(ax.get_position().bounds)
        x_text = pos[0] - 0.01
        y_text = pos[1] + pos[3] / 2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axes:
        ax.set_axis_off()

    plt.savefig('colormaps.png', dpi=100)
