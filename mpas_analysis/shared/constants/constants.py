import numpy as np

"""
	Constants that are common to all ocean model vs observations analysis

	Luke Van Roekel
	10/21/2016

"""

#set parameters for interpolated grid
dLongitude = 1.
dLatitude = 1.
lonmin = -180.
lonmax = 180.
latmin = -90.
latmax = 90.

monthdictionary={'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6, 'Jul':7, 'Aug':8, 'Sep':9, 'Oct':10,
                 'Nov':11, 'Dec':12, 'JFM':np.array([1,2,3]), 'AMJ':np.array([4,5,6]), 'JAS':np.array([7,8,9]),
                 'OND':np.array([10,11,12]), 'ANN':np.arange(1,13)}

dinmonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])


