# make sure to set the Agg backend that works even without x-forwarding
# before any other matplotlib imports
import sys
import warnings
import matplotlib as mpl
mpl.use('Agg')

if sys.version_info[0] < 3:
    print("Warning: MPAS-Analysis v1.2.5 is the last to support python 2.7")

__version_info__ = (1, 2, 5)
__version__ = '.'.join(str(vi) for vi in __version_info__)
