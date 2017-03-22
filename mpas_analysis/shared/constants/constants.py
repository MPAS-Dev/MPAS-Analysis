import numpy as np

"""
Constants that are common to all analysis tasks

Authors
-------
Luke Van Roekel, Xylar Asay-Davis, Milena Veneziani

Last modified
-------------
03/15/2017
"""

# set parameters for default climatology comparison grid
dLongitude = 0.5
dLatitude = 0.5
lonmin = -180.
lonmax = 180.
latmin = -90.
latmax = 90.

monthDictionary = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
                   'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11,
                   'Dec': 12, 'JFM': np.array([1, 2, 3]),
                   'AMJ': np.array([4, 5, 6]), 'JAS': np.array([7, 8, 9]),
                   'OND': np.array([10, 11, 12]), 'ANN': np.arange(1, 13),
                   'ON': np.array([10, 11]), 'FM': np.array([2, 3]),
                   'DJF': np.array([12, 1, 2]), 'JJA': np.array([6, 7, 8])}

daysInMonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

abrevMonthNames = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]

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

# small value to prevent division by zero
eps = 1.E-10

# vim: foldmethod=marker ai ts=4 sts=4 et sw=4 ft=python
