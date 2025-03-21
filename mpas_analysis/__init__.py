# make sure to set the Agg backend that works even without x-forwarding
# before any other matplotlib imports
import matplotlib as mpl
mpl.use('Agg')
