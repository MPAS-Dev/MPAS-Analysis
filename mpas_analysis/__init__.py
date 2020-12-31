# make sure to set the Agg backend that works even without x-forwarding
# before any other matplotlib imports
import matplotlib as mpl
mpl.use('Agg')

__version_info__ = (1, 3, 0)
__version__ = '.'.join(str(vi) for vi in __version_info__)
