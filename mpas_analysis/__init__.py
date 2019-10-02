# make sure to set the Agg backend that works even without x-forwarding
# before any other matplotlib imports
import sys
import warnings
import matplotlib as mpl
mpl.use('Agg')

__version_info__ = (1, 2, 5)
__version__ = '.'.join(str(vi) for vi in __version_info__)
