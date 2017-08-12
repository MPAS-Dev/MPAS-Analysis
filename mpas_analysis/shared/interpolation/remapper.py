'''
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
04/13/2017
'''

import subprocess
import tempfile
import os
from distutils.spawn import find_executable
import numpy
from scipy.sparse import csr_matrix
import xarray as xr
import sys

from ..grid import MpasMeshDescriptor, LatLonGridDescriptor, \
    ProjectionGridDescriptor


class Remapper(object):
    '''
    A class for remapping fields using a given mapping file.  The weights and
    indices from the mapping file can be loaded once and reused multiple times
    to map several fields between the same source and destination grids.
    Author
    ------
    Xylar Asay-Davis

    Last Modified
    -------------
    04/13/2017
    '''

    def __init__(self, sourceDescriptor, destinationDescriptor,
                 mappingFileName=None):  # {{{
        '''
        Create the remapper and read weights and indices from the given file
        for later used in remapping fields.

        Parameters
        ----------
        sourceDescriptor : an instance of {MpasMeshDescriptor,
                                           LatLonGridDescriptor,
                                           ProjectionGridDescriptor}
            An object used to write a scrip file and to determine the type of
            the source mesh or grid.

        destinationDescriptor : an instance of {MpasMeshDescriptor,
                                                LatLonGridDescriptor,
                                                ProjectionGridDescriptor}
            An object used to write a scrip files and to determine the type of
            the destination mesh or grid.

        mappingFileName : str, optional
            The path where the mapping file containing interpolation weights
            and indices will be written and/or read.  If ``None``,
            no interpolation is performed and data sets are returned unchanged.
            This is useful if the source and destination grids are determined
            to be the same (though the Remapper does not attempt to determine
            if this is the case).

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/13/2017
        '''

        if not isinstance(sourceDescriptor,
                          (MpasMeshDescriptor,  LatLonGridDescriptor,
                           ProjectionGridDescriptor)):
            raise ValueError("sourceDescriptor is not of a recognized type.")

        if not isinstance(destinationDescriptor,
                          (MpasMeshDescriptor,  LatLonGridDescriptor,
                           ProjectionGridDescriptor)):
            raise ValueError(
                "destinationDescriptor is not of a recognized type.")

        self.sourceDescriptor = sourceDescriptor
        self.destinationDescriptor = destinationDescriptor
        self.mappingFileName = mappingFileName

        self.mappingLoaded = False

        # }}}

    def build_mapping_file(self, method='bilinear',
                           additionalArgs=None):  # {{{
        '''
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        method : {'bilinear', 'neareststod', 'conserve'}, optional
            The method of interpolation used, see documentation for
            `ESMF_RegridWeightGen` for details.

        additionalArgs : list of str, optional
            A list of additional arguments to ``ESMF_RegridWeightGen``

        Raises
        ------
        OSError
            If ``ESMF_RegridWeightGen`` is not in the system path.

        ValueError
            If sourceDescriptor or destinationDescriptor is of an unknown type

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/13/2017
        '''

        if self.mappingFileName is None or \
                os.path.exists(self.mappingFileName):
            # a valid weight file already exists, so nothing to do
            return

        if find_executable('ESMF_RegridWeightGen') is None:
            raise OSError('ESMF_RegridWeightGen not found. Make sure esmf '
                          'package is installed via\n'
                          'latest nco: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        # Write source and destination SCRIP files in temporary locations
        self.sourceDescriptor.to_scrip(_get_temp_path())
        self.destinationDescriptor.to_scrip(_get_temp_path())

        args = ['ESMF_RegridWeightGen',
                '--source', self.sourceDescriptor.scripFileName,
                '--destination', self.destinationDescriptor.scripFileName,
                '--weight', self.mappingFileName,
                '--method', method,
                '--netcdf4',
                '--no_log']

        if self.sourceDescriptor.regional:
            args.append('--src_regional')

        if self.destinationDescriptor.regional:
            args.append('--dst_regional')

        if self.sourceDescriptor.regional or \
                self.destinationDescriptor.regional:
            args.append('--ignore_unmapped')

        if additionalArgs is not None:
            args.extend(additionalArgs)

        # make sure any output is flushed before we add output from the
        # subprocess
        sys.stdout.flush()
        sys.stderr.flush()
        subprocess.check_call(args)

        # remove the temporary SCRIP files
        os.remove(self.sourceDescriptor.scripFileName)
        os.remove(self.destinationDescriptor.scripFileName)

        # }}}

    def remap_file(self, inFileName, outFileName,
                   variableList=None, overwrite=False):  # {{{
        '''
        Given a source file defining either an MPAS mesh or a lat-lon grid and
        a destination file or set of arrays defining a lat-lon grid, constructs
        a mapping file used for interpolation between the source and
        destination grids.

        Parameters
        ----------
        inFileName : str
            The path to the file containing a data set on the source grid

        outFileName : str
            The path where the data on the destination grid should be written

        variableList : list of str, optional
            A list of variables to be mapped.  By default, all variables are
            mapped

        overwrite : bool, optional
            Whether the destination file should be overwritten if it already
            exists. If `False`, and the destination file is already present,
            the function does nothing and returns immediately

        Raises
        ------
        OSError
            If ``ncremap`` is not in the system path.

        ValueError
            If ``mappingFileName`` is ``None`` (meaning no remapping is
            needed).

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/13/2017
        '''

        if self.mappingFileName is None:
            raise ValueError('No mapping file was given because remapping is '
                             'not necessary. The calling\n'
                             'code should simply use the constents of {} '
                             'directly.'.format(inFileName))

        if not overwrite and os.path.exists(outFileName):
            # a remapped file already exists, so nothing to do
            return

        if isinstance(self.sourceDescriptor, ProjectionGridDescriptor):
            raise TypeError('Source grid is a projection grid, not supported '
                            'by ncremap.\n'
                            'Consider using Remapper.remap')
        if isinstance(self.destinationDescriptor, ProjectionGridDescriptor):
            raise TypeError('Destination grid is a projection grid, not '
                            'supported by ncremap.\n'
                            'Consider using Remapper.remap')

        if find_executable('ncremap') is None:
            raise OSError('ncremap not found. Make sure the latest nco '
                          'package is installed: \n'
                          'conda install nco\n'
                          'Note: this presumes use of the conda-forge '
                          'channel.')

        args = ['ncremap',
                '-i', inFileName,
                '-m', self.mappingFileName,
                '--vrb=1',
                '-o', outFileName]

        if isinstance(self.sourceDescriptor, LatLonGridDescriptor):
            args.extend(['-R', '--rgr lat_nm={}  --rgr lon_nm={}'.format(
                self.sourceDescriptor.latVarName,
                self.sourceDescriptor.lonVarName)])

        if isinstance(self.sourceDescriptor, MpasMeshDescriptor):
            # Note: using the -C (climatology) flag for now because otherwise
            #       ncremap tries to add a _FillValue attribute that might
            #       already be present and quits with an error
            args.extend(['-P', 'mpas', '-C'])

        if variableList is not None:
            args.extend(['-v', ','.join(variableList)])

        # make sure any output is flushed before we add output from the
        # subprocess
        sys.stdout.flush()
        sys.stderr.flush()

        subprocess.check_call(args)  # }}}

    def remap(self, ds, renormalizationThreshold=None):  # {{{
        '''
        Given a source data set, returns a remapped version of the data set,
        possibly masked and renormalized.

        Parameters
        ----------
        ds : ``xarray.Dataset`` or ``xarray.DataArray`` object
            The dimention(s) along ``self.sourceDimNames`` must match
            ``self.src_grid_dims`` read from the mapping file.

        renormalizationThreshold : float, optional
            The minimum weight of a denstination cell after remapping, below
            which it is masked out, or ``None`` for no renormalization and
            masking.

        Returns
        -------
        remappedDs : `xarray.Dataset`` or ``xarray.DataArray`` object
            Returns a remapped data set (or data array) where dimensions other
            than ``self.sourceDimNames`` are the same as in ``ds`` and the
            dimension(s) given by ``self.sourceDimNames`` have been replaced by
            ``self.destinationDimNames``.

        Raises
        ------
        ValueError
            If the size of ``self.sourceDimNames`` in ``ds`` do not match the
            source dimensions read in from the mapping file
            (``self.src_grid_dims``).
        TypeError
            If ds is not an ``xarray.Dataset`` or ``xarray.DataArray`` object

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/13/2017
        '''

        if self.mappingFileName is None:
            # No remapping is needed
            return ds

        self._load_mapping()

        for index, dim in enumerate(self.sourceDescriptor.dims):
            if self.src_grid_dims[index] != ds.sizes[dim]:
                raise ValueError('data set and remapping source dimension {} '
                                 'don\'t have the same size: {} != {}'.format(
                                     dim, self.src_grid_dims[index],
                                     len(ds.sizes[dim])))

        if isinstance(ds, xr.DataArray):
            remappedDs = self._remap_data_array(ds, renormalizationThreshold)
        elif isinstance(ds, xr.Dataset):
            drop = []
            for var in ds.data_vars:
                if self._check_drop(ds[var]):
                    drop.append(var)
            remappedDs = ds.drop(drop)
            remappedDs = remappedDs.apply(self._remap_data_array,
                                          keep_attrs=True,
                                          args=(renormalizationThreshold,))
        else:
            raise TypeError('ds not an xarray Dataset or DataArray.')

        # Update history attribute of netCDF file
        if 'history' in remappedDs.attrs:
            newhist = '\n'.join([remappedDs.attrs['history'],
                                 ' '.join(sys.argv[:])])
        else:
            newhist = sys.argv[:]
        remappedDs.attrs['history'] = newhist

        remappedDs.attrs['meshName'] = self.destinationDescriptor.meshName

        return remappedDs  # }}}

    def _load_mapping(self):  # {{{
        '''
        Load weights and indices from a mapping file, if this has not already
        been done

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/06/2017
        '''

        if self.mappingLoaded:
            return

        dsMapping = xr.open_dataset(self.mappingFileName)
        n_a = dsMapping.dims['n_a']
        n_b = dsMapping.dims['n_b']

        nSourceDims = len(self.sourceDescriptor.dims)
        src_grid_rank = dsMapping.dims['src_grid_rank']
        nDestinationDims = len(self.destinationDescriptor.dims)
        dst_grid_rank = dsMapping.dims['dst_grid_rank']

        # check that the mapping file has the right number of dimensions
        if nSourceDims != src_grid_rank or \
                nDestinationDims != dst_grid_rank:
            raise ValueError('The number of source and/or '
                             'destination dimensions does not\n'
                             'match the expected number of source and '
                             'destination dimensions in the mapping\n'
                             'file. {} != {} and/or {} != {}'.format(
                                 nSourceDims, src_grid_rank,
                                 nDestinationDims, dst_grid_rank))

        # grid dimensions need to be reversed because they are in Fortran order
        self.src_grid_dims = dsMapping['src_grid_dims'].values[::-1]
        self.dst_grid_dims = dsMapping['dst_grid_dims'].values[::-1]

        # now, check that each source and destination dimension is right
        for index in range(len(self.sourceDescriptor.dims)):
            dim = self.sourceDescriptor.dims[index]
            dimSize = self.sourceDescriptor.dimSize[index]
            checkDimSize = self.src_grid_dims[index]
            if dimSize != checkDimSize:
                raise ValueError('source mesh descriptor and remapping source '
                                 'dimension {} don\'t have the same size: \n'
                                 '{} != {}'.format(dim, dimSize, checkDimSize))
        for index in range(len(self.destinationDescriptor.dims)):
            dim = self.destinationDescriptor.dims[index]
            dimSize = self.destinationDescriptor.dimSize[index]
            checkDimSize = self.dst_grid_dims[index]
            if dimSize != checkDimSize:
                raise ValueError('dest. mesh descriptor and remapping dest. '
                                 'dimension {} don\'t have the same size: \n'
                                 '{} != {}'.format(dim, dimSize, checkDimSize))

        self.frac_b = dsMapping['frac_b'].values

        col = dsMapping['col'].values-1
        row = dsMapping['row'].values-1
        S = dsMapping['S'].values
        self.matrix = csr_matrix((S, (row, col)), shape=(n_b, n_a))

        self.mappingLoaded = True  # }}}

    def _check_drop(self, dataArray):  # {{{
        sourceDims = self.sourceDescriptor.dims

        sourceDimsInArray = [dim in dataArray.dims for dim in sourceDims]

        return (numpy.any(sourceDimsInArray) and not
                numpy.all(sourceDimsInArray))  # }}}

    def _remap_data_array(self, dataArray, renormalizationThreshold):  # {{{
        '''
        Regrids a single xarray data array

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/05/2017
        '''

        sourceDims = self.sourceDescriptor.dims
        destDims = self.destinationDescriptor.dims

        sourceDimsInArray = [dim in dataArray.dims for dim in sourceDims]

        if not numpy.any(sourceDimsInArray):
            # no remapping is needed
            return dataArray

        if not numpy.all(sourceDimsInArray):
            # no remapping is possible so the variable array should have been
            # dropped
            raise ValueError('Data array with some (but not all) required '
                             'source dims cannot be remapped\n'
                             'and should have been dropped.')

        # make a list of dims and remapAxes
        dims = []
        remapAxes = []
        destDimsAdded = False
        for index, dim in enumerate(dataArray.dims):
            if dim in sourceDims:
                remapAxes.append(index)
                if not destDimsAdded:
                    dims.extend(destDims)
                    destDimsAdded = True
            else:
                dims.append(dim)

        # make a dict of coords
        coordDict = {}
        # copy unmodified coords
        for coord in dataArray.coords:
            sourceDimInCoord = numpy.any([dim in dataArray.coords[coord].dims
                                          for dim in sourceDims])
            if not sourceDimInCoord:
                coordDict[coord] = {'dims': dataArray.coords[coord].dims,
                                    'data': dataArray.coords[coord].values}

        # add dest coords
        coordDict.update(self.destinationDescriptor.coords)

        # remap the values
        field = dataArray.values
        mask = numpy.isnan(field)
        if numpy.count_nonzero(mask) > 0:
            field = numpy.ma.masked_array(field, mask)
        remappedField = self._remap_numpy_array(field, remapAxes,
                                                renormalizationThreshold)

        arrayDict = {'coords': coordDict,
                     'attrs': dataArray.attrs,
                     'dims': dims,
                     'data': remappedField,
                     'name': dataArray.name}

        # make a new data array
        remappedArray = xr.DataArray.from_dict(arrayDict)

        return remappedArray  # }}}

    def _remap_numpy_array(self, inField, remapAxes,
                           renormalizationThreshold):  # {{{
        '''
        Regrids a single numpy array

        Author
        ------
        Xylar Asay-Davis

        Last Modified
        -------------
        04/05/2017
        '''

        # permute the dimensions of inField so the axes to remap are first,
        # then flatten the remapping and the extra dimensions separately for
        # the matrix multiply
        extraAxes = [axis for axis in numpy.arange(inField.ndim)
                     if axis not in remapAxes]

        newShape = [numpy.prod([inField.shape[axis] for axis in remapAxes])]
        if len(extraAxes) > 0:
            extraShape = [inField.shape[axis] for axis in extraAxes]
            newShape.append(numpy.prod(extraShape))
        else:
            extraShape = []
            newShape.append(1)

        permutedAxes = remapAxes + extraAxes

        # permute axes so the remapped dimension(s) come first and "flatten"
        # the remapping dimension
        inField = inField.transpose(permutedAxes).reshape(newShape)

        masked = (isinstance(inField, numpy.ma.MaskedArray) and
                  renormalizationThreshold is not None)
        if masked:
            inMask = numpy.array(numpy.logical_not(inField.mask), float)
            outField = self.matrix.dot(inMask*inField)
            outMask = self.matrix.dot(inMask)
            mask = outMask > renormalizationThreshold
        else:
            outField = self.matrix.dot(inField)
            # make frac_b match the shape of outField
            outMask = numpy.reshape(self.frac_b, (len(self.frac_b), 1)).repeat(
                newShape[1], axis=1)
            mask = outMask > 0.

        # normalize the result based on outMask
        outField[mask] /= outMask[mask]
        outField = numpy.ma.masked_array(outField,
                                         mask=numpy.logical_not(mask))

        destRemapDimCount = len(self.dst_grid_dims)
        outDimCount = len(extraShape) + destRemapDimCount

        # "unflatten" the remapped dimension(s)
        destShape = list(self.dst_grid_dims) + extraShape
        outField = numpy.reshape(outField, destShape)

        # "unpermute" the axes to be in the expected order
        index = numpy.amin(remapAxes)
        unpermuteAxes = list(numpy.arange(destRemapDimCount, outDimCount))
        unpermuteAxes = (unpermuteAxes[0:index] +
                         list(numpy.arange(destRemapDimCount)) +
                         unpermuteAxes[index:])
        outField = numpy.transpose(outField, axes=unpermuteAxes)

        return outField  # }}}


def _get_lock_path(fileName):  # {{{
    '''Returns the name of a temporary lock file unique to a given file name'''
    directory = '{}/.locks/'.format(os.path.dirname(fileName))
    try:
        os.makedirs(directory)
    except OSError:
        pass
    return '{}/{}.lock'.format(directory, os.path.basename(fileName))  # }}}


def _get_temp_path():  # {{{
    '''Returns the name of a temporary NetCDF file'''
    return '{}/{}.nc'.format(tempfile._get_default_tempdir(),
                             next(tempfile._get_candidate_names()))  # }}}

# vim: ai ts=4 sts=4 et sw=4 ft=python
