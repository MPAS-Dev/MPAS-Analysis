"""
    ScientificColourMaps5

    Usage
    -----
    import ScientificColourMaps5 as SCM5
    plt.imshow(data, cmap=SCM5.berlin)

    Available colourmaps
    ---------------------
    acton, bamako, batlow, berlin, bilbao, broc, buda, cork, davos, devon,
    grayC, hawaii, imola, lajolla, lapaz, lisbon, nuuk, oleron, oslo, roma,
    tofino, tokyo, turku, vik
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

folder = os.path.abspath(os.path.dirname(os.path.abspath(__file__)))

__all__ = {'acton', 'bamako', 'batlow', 'berlin', 'bilbao', 'broc', 'buda',
           'cork', 'davos', 'devon', 'grayC', 'hawaii', 'imola', 'lajolla',
           'lapaz', 'lisbon', 'nuuk', 'oleron', 'oslo', 'roma', 'tofino',
           'tokyo', 'turku', 'vik'}

for name in __all__:
    file = os.path.join(folder, name, name + '.txt')
    cm_data = np.loadtxt(file)
    cmap = LinearSegmentedColormap.from_list(name, cm_data)
    vars()[name] = cmap
    plt.register_cmap(name, cmap)
    plt.register_cmap('{}_r'.format(name), cmap.reversed())