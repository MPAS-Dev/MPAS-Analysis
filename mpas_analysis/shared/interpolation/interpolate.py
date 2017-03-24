"""
Functions for performing interpolation

Functions
---------
build_remap_weights - constructs a mapping file containing the indices and
    weights needed to perform horizontal interpolation

remap - perform horizontal interpolation on a data sets, given a mapping file

Author
------
Xylar Asay-Davis

Last Modified
-------------
03/14/2017
"""

import subprocess
import tempfile
import os
from distutils.spawn import find_executable

from .scrip import mpas_file_to_scrip, lat_lon_file_to_scrip, \
    lat_lon_array_to_scrip


def build_remap_weights(sourceFileName, outWeightFileName,
                        destintionFileName=None, sourceFileType='mpas',
                        sourceLatVarName='lat', sourceLonVarName='lon',
                        destintionLatVarName='lat', destintionLonVarName='lon',
                        destinationLat=None, destinationLon=None,
                        desitnationUnits='degrees',
                        method='bilinear', overwrite=False):  # {{{
    """
    Given a source file defining either an MPAS mesh or a lat-lon grid and
    a destination file or set of arrays defining a lat-lon grid, constructs
    a mapping file used for interpolation between the source and destination
    grids.

    Parameters
    ----------
    sourceFileName : str
        The path of the file containing either the source MPAS mesh or
        the source lat-lon grid

    outWeightFileName : str
        The path to which the mapping file containing interpolation weights
        and indices should be written

    destintionFileName : str, optional
        The path of the file containing the destination lat-lon grid.  Should
        be None if `destinationLat` and `destinationLon` are supplied instead.

    sourceFileType : {'mpas', 'latlon'}
        Whether the source file contains an MPAS mesh or a lat-lon grid

    sourceLatVarName, sourceLonVarName : str, optional
        If `sourceFileType == 'latlon'`, the name of the latitude and longitude
        variables in the source grid file

    destintionLatVarName, destintionLonVarName : str, optional
        If `destintionFileName` is not `None`, the name of the latitude and
        longitude variables in the source grid file

    destinationLat, destinationLon : 1D numpy.arrays, optional
        One dimensional arrays defining the latitude and longitude coordinates
        of grid corners on the destination grid. `destintionFileName` should be
        set to `None` if these are supplied

    desitnationUnits : {'degrees', 'radians'}, optional
        The units of `destinationLat` and `destinationLon` (if they are
        supplied)

    method : {'bilinear', 'neareststod', 'conserve'}
        The method of interpolation used, see documentation for
        `ESMF_RegridWeightGen` for details.

    overwrite : bool, optional
        Whether the mapping file should be overwritten if it already exists.
        If `False`, and the mapping file is already present, the function
        does nothing and returns immediately, potentially saving a costly
        re-computaiton of the mapping file.

    Raises
    ------
    OSError
        If `ESMF_RegridWeightGen` is not in the system path.

    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/14/2017
    """

    if not overwrite and os.path.exists(outWeightFileName):
        # a valid weight file already exists, so nothing to do
        return

    if find_executable('ESMF_RegridWeightGen') is None:
        raise OSError('ESMF_RegridWeightGen not found. Make sure esmf package '
                      'is installed via\nlatest nco: \n'
                      'conda install nco\n'
                      'Note: this presumes use of the conda-forge channel.')

    # two temporary SCRIP files, one for the MPAS mesh and one for the dest
    # grid
    sourceScripFileName = _get_temp_path()
    destintionScripFileName = _get_temp_path()

    args = ['ESMF_RegridWeightGen', '--source', sourceScripFileName,
            '--destination', destintionScripFileName,
            '--weight', outWeightFileName,
            '--method', method]

    if sourceFileType == 'mpas':
        mpas_file_to_scrip(mpasFileName=sourceFileName,
                           scripFileName=sourceScripFileName)
        args.extend(['--src_regional', '--ignore_unmapped'])
    elif sourceFileType == 'latlon':
        lat_lon_file_to_scrip(inFileName=sourceFileName,
                              scripFileName=sourceScripFileName,
                              latVarName=sourceLatVarName,
                              lonVarName=sourceLonVarName)
    else:
        raise ValueError("sourceFileType is neither 'mpas' or 'latlon'.")

    if destintionFileName is not None:
        lat_lon_file_to_scrip(inFileName=destintionFileName,
                              scripFileName=destintionScripFileName,
                              latVarName=destintionLatVarName,
                              lonVarName=destintionLonVarName)
    elif destinationLat is not None and destinationLon is not None:
        lat_lon_array_to_scrip(latCorner=destinationLat,
                               lonCorner=destinationLon,
                               units=desitnationUnits,
                               scripFileName=destintionScripFileName)
    else:
        raise ValueError('Either destintionFileName or both config and '
                         'sectionName must be supplied.')

    subprocess.check_call(args)

    # remove the temporary SCRIP files
    os.remove(sourceScripFileName)
    os.remove(destintionScripFileName)  # }}}


def remap(inFileName, outFileName, inWeightFileName, sourceFileType='mpas',
          sourceLatVarName='lat', sourceLonVarName='lon',
          variableList=None, overwrite=False):  # {{{
    """
    Given a source file defining either an MPAS mesh or a lat-lon grid and
    a destination file or set of arrays defining a lat-lon grid, constructs
    a mapping file used for interpolation between the source and destination
    grids.

    Parameters
    ----------
    inFileName : str
        The path to the file containing a data set on the source grid

    outFileName : str
        The path where the data on the destination grid should be written

    inWeightFileName : str
        The path to the mapping file containing interpolation weights
        and indices between the source and destination grids

    sourceFileType : {'mpas', 'latlon'}
        Whether the source file contains an MPAS mesh or a lat-lon grid

    sourceLatVarName, sourceLonVarName : str, optional
        If `sourceFileType == 'latlon'`, the name of the latitude and longitude
        variables in the source grid file

    variableList : list of str, optional
        A list of variables to be mapped.  By default, all variables are mapped

    overwrite : bool, optional
        Whether the destination file should be overwritten if it already
        exists. If `False`, and the destination file is already present, the
        function does nothing and returns immediately

    Raises
    ------
    OSError
        If `ncremap` is not in the system path.

    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    03/14/2017
    """

    if not overwrite and os.path.exists(outFileName):
        # a valid weight file already exists, so nothing to do
        return

    if find_executable('ncremap') is None:
        raise OSError('ncremap not found. Make sure the latest nco package '
                      'is installed: \n    conda install nco')

    args = ['ncremap',
            '-R', '--rgr lat_nm={} --rgr lon_nm={}'.format(sourceLatVarName,
                                                           sourceLonVarName),
            '-i', inFileName,
            '-m', inWeightFileName,
            '-o', outFileName]

    if sourceFileType == 'mpas':
        # Note: using the -C (climatology) flag for now because otherwise
        #       ncremap tries to add a _FillValue attribute that might already
        #       be present and quits with an error
        args.extend(['-P', 'mpas', '-C'])
    if variableList is not None:
        args.extend(['-v', ','.join(variableList)])

    subprocess.check_call(args)  # }}}


def _get_temp_path():  # {{{
    '''Returns the name of a temporary NetCDF file'''
    return '{}/{}.nc'.format(tempfile._get_default_tempdir(),
                             next(tempfile._get_candidate_names()))  # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
