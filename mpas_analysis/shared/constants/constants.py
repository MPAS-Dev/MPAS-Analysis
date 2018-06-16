"""
Constants that are common to all analysis tasks
"""
# Authors
# -------
# Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

from __future__ import absolute_import, division, print_function, \
    unicode_literals

import numpy as np

# set parameters for default climatology comparison grid
dLongitude = 0.5
dLatitude = 0.5
lonmin = -180.
lonmax = 180.
latmin = -90.
latmax = 90.

monthsInYear = 12

monthDictionary = {'Jan': [1], 'Feb': [2], 'Mar': [3], 'Apr': [4], 'May': [5],
                   'Jun': [6], 'Jul': [7], 'Aug': [8], 'Sep': [9], 'Oct': [10],
                   'Nov': [11], 'Dec': [12], 'JFM': [1, 2, 3],
                   'AMJ': [4, 5, 6], 'JAS': [7, 8, 9], 'OND': [10, 11, 12],
                   'ANN': list(np.arange(1, 13)), 'ON': [10, 11], 'FM': [2, 3],
                   'DJF': [12, 1, 2], 'JJA': [6, 7, 8]}

daysInMonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

abrevMonthNames = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug",
                   "Sep", "Oct", "Nov", "Dec"]

# conversion factor from m^3/s to Sverdrups
m3ps_to_Sv = 1e-6

# conversion factor from radians to degrees
rad_to_deg = 180./np.pi

# conversion factor from degrees to radians
deg_to_rad = np.pi/180.

# seconds in a year
sec_per_year = 86400. * 365.

# seconds per month (approximate)
sec_per_month = 86400. * 30.

# Tapering coefficient for calculating spectral degrees of freedom
tapcoef = 1.055111111111111

# small value to prevent division by zero
eps = 1.E-10

# density of freshwater (kg/m^3)
rho_fw = 1000.

# kilograms per gigatonne
kg_per_GT = 1e12

# cm per m
cm_per_m = 100.

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
