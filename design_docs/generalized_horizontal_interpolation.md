<h1> Title: Generalized Horizontal Interpolation in MPAS-Analysis <br>
Xylar Asay-Davis <br>
date: 01-19-2017 <br>
</h1>
<h2> Summary </h2>
Currently, MPAS-Analysis uses various methods to perform horizontal interpolation.  For constructing ocean climatologies, nearest-neighbor interpolation is used, while for sea-ice climatologies, `ncremap` is used with the requirement that a mapping file for the appropriate source and destination grids is provided through the config file.  This project intends to move MPAS-Analysis to a unified approach to horizontal interpolation that does not require pre-generated mapping files (though it should support caching mapping files for faster execution).

Many types of analysis in MPAS will require fields that are interpolated from MPAS grids to arbitrary points, not just to points on a lat/lon grid.  This project will not attempt to address that case completely but will take that need into consideration in designing a solution that can be extended to interpolation at arbitrary points in the future.

<h1> Requirements </h1>
<h2>Requirement: Higher-order interpolation <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

The option to interpolate smoothly (e.g. linearly or with barycentric coordinates) between cell-centered values should be added.  The calling code should easily be able to select among various orders of interpolation with a flag.

<h2>Consideration: Interpolation should handle periodic boundaries <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

If and when MPAS-Analysis supports planar test cases with periodic boundaries, interpolation should be extended to handle periodic boundaries

<h2>Consideration: Interpolation should handle Cartesian meshes <br>
Date last modified:  2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

If and when MPAS-Analysis supports planar test cases with purely Cartesian meshes (e.g. where `latCell` and `lonCell` do not vary), interpolation should be extended to handle Cartesian Coordinates

<h2>Consideration: Support for arbitrary output interpolation points <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis
</h2>

The calling code should be able to supply any desired interpolation points, not just a regular latitude-longitude grid.

<h2>Consideration: Support caching results from any costly, one-time geometric computations <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

For many potential algorithms used to perform interpolation, there is likely to be a relatively costly step of computing fields such as indices into input data fields and interpolation weights that 1) only need to be computed once for a given input mesh and set of output points and 2) are independent of the data in the field being interpolated.  If this data were cached, it could mean that rerunning the analysis (which might be very desirable, e.g., while monitoring the progress of a run) would be much cheaper than the initial run.  Also, a cached weight file from a previous analysis run could be used when analyzing a subsequent run with identical source meshes.



<h1> Algorithmic Formulations </h1>

<h2>Design solution: Higher-order interpolation <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

The approach will be to create SCRIP files (or, in the future for greater flexibility perhaps ESMF grid/mesh files) for the source and destination grids, then to use `ESMF_RegridWeightGen` to generate a mapping file.  `ESMF_RegridWeightGen` supports 5 interpolation methods---bilinear, patch, nearestdtos, neareststod, and conserve---and we would likely support at least bilinear, neareststod and conserve, and perhaps all 5.  The destination grid will be specified either by reading values from `lat` and `lon` coordinates of a NetCDF file or through config file options `lat` and `lon` that are typically expressions involving `numpy.arange` or `numpy.linspace`.

Then, `ncremap` will be used to remap the desired list of variables from an MPAS NetCDF file to the desired destination grid.

<h2>Design solution: Interpolation should handle periodic boundaries <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

For now, periodic boundaries (except for the obvious one at +/- 180 longitude) will not be supported.  It appears that ESMF grid files do include support for periodic boundaries so the current solution should be relatively easy to extend to periodic boundaries in the future.

<h2>Design solution: Interpolation should handle Cartesian meshes <br>
Date last modified:  2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

ESMF unstructured mesh files seem to support Cartesian coordinates.  This will be investigated if and when MPAS-Analysis can accommodate a test case with Cartesian coordinates.

<h2>Design solution: Support for arbitrary output interpolation points <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis
</h2>

I do not intend to address this consideration in this project.  It may be that `ESMF_RegridWeightGen` can also be used to perform interpolation to arbitrary points (in particular, a set of points that are not cell centers or vertices of a mesh), but this is not yet clear to me.  If not, an alternative solution for arbitrary destination points will be needed.

<h2>Design solution: Support caching results from any costly, one-time geometric computations <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

This should be relatively easy to accommodate with `ESMF_RegridWeightGen` and `ncremap`.  The default behavior of the function for generating interpolation weights will be to do nothing if the mapping file already exists.  Further, we can support an optional config option that will point to an existing mapping file if one has already been generated and cached somewhere (e.g. in a shared directory).  Eventually, we will probably want to systematically store these mapping files for typical MPAS meshes and typical output grids, particularly for those that are expensive to generate.

<h1> Design and Implementation </h1>
<h2>Implementation: Higher-order interpolation <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Implementation is in the branch https://github.com/xylar/MPAS-Analysis/tree/horiz_interp.

`ESMF_RegridWeightGen` is used to compute regridding weights that are 'bilinear', 'neareststod' (nearest neighbor) or 'conserve' (conservative).  The order of regridding can be chosen separately for MPAS model results, ocean observations and sea-ice observationos via `mpasInterpolationMethod` and `interpolationMethod` flags (see the template: https://github.com/xylar/MPAS-Analysis/blob/horiz_interp/config.template).

<h2>Implementation: Interpolation should handle periodic boundaries <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Not yet supported.

<h2>Implementation: Interpolation should handle Cartesian meshes <br>
Date last modified:  2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Not yet supported.

<h2>Implementation: Support for arbitrary output interpolation points <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis
</h2>

Not yet supported.

<h2>Implementation: Support caching results from any costly, one-time geometric computations <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

Mapping files, climatologies and remapped climatologies are cached when they are created.  Both mapping files and the directory containing the remapped climatologies from observations can be supplied via the config file, saving the time of computing them.



<h1> Testing </h1>
<h2>Testing and Validation: Higher-order interpolation <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Testing of each of the flags ('bilinear', 'neareststod' and 'conserve') has been performed with the `GMPAS-QU240` run, all of wich produce plots that look acceptable.  Bilinear and conserve methods leave halos of invalid cells around land at coarse resolution, which is consistent with the coarse resolution of this test mesh.

An alpha8 and a beta0 run on edison.  They ran successfully but I have not had a chance to examine the output.

<h2>Testing and Validation: Interpolation should handle periodic boundaries <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Not yet supported.

<h2>Testing and Validation: Interpolation should handle Cartesian meshes <br>
Date last modified:  2017/03/04 <br>
Contributors: Xylar Asay-Davis</h2>

Not yet supported.

<h2>Testing and Validation: Support for arbitrary output interpolation points <br>
Date last modified: 2017/03/04 <br>
Contributors: Xylar Asay-Davis
</h2>

Not yet supported.

<h2>Testing and Validation: Support caching results from any costly, one-time geometric computations <br>
Date last modified: 2017/02/25 <br>
Contributors: Xylar Asay-Davis</h2>

I have verified that I can rerun without re-computing mapping files or climatologies.  Using the `GMPAS-QU240` run, I have verified that I can supply mapping files and remapped observation climatologies without them being re-computed
