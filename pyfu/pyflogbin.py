# Copyright(c) 2013-2014 Association of Universities for Research in Astronomy, Inc.,
# by James E.H. Turner.
#
# 'pyflogbin' Python module for rebinning cubes logarithmically in wavelength
#
# Version  Nov, 2013  JT Initial version

import pyfits, numpy
import scipy.ndimage.interpolation as ndi
import astro_ds

import cProfile, pstats, sys

reload(astro_ds)

# Pyflogbin main routine (non-PyRAF interface):
def pyflogbin(inimage, outimage, fluxcons=False, propvar=False):

    """Re-bin a data cube with a logarithmic spectral WCS"""

    # Open the input file:
    inhdulist = pyfits.open(inimage, mode='readonly')
    inds = astro_ds.DataSet(inhdulist, 1)   # Just assume SCI,1 for now

    # Copy the output PHU from the input file:
    outphu = pyfits.PrimaryHDU(header=inhdulist[0].header.copy())

    # Save the PHU to disk, to ensure the file is writeable:
    outphu.writeto(outimage)
    del outphu

    # Re-open the minimal output image for appending:
    outhdulist = pyfits.open(outimage, mode='append')

    # Create an empty DataSet for the output and set the WCS and array
    # size, based on the WCS and array sizes of the list of input datasets:
    outds = astro_ds.DataSet()
    outds.SetGridFrom([inds], wavtolog=True)

    # Specify log wavelength sampling:
    ReBinWavelength(inds, outds, fluxcons=False, propvar=propvar)

    # Copy the SCI header from the input file as a starting point
    # for the output header:
    scihdr = inhdulist['sci',1].header.copy()
    try:
        varhdr = inhdulist['var',1].header.copy()
    except KeyError:
        varhdr = None

    # Append the output DataSet as an extension of the output HDUList:
    outds.WriteHDU(outhdulist, header=scihdr, varheader=varhdr)

    # Close files, writing the output to disk:
    outhdulist.close(output_verify="warn")
    inhdulist.close(output_verify="warn")

# End (pyflogbin main Python routine)


def ReBinWavelength(inds, outds, fluxcons=False, propvar=False):
    """Resample a DataSet with a previously-defined output wavelength binning"""

    # Propagate variance if asked to and it's available in the input:
    if propvar:
        propvar = inds.GetVar(create=False) is not None

    incube = inds.GetData()
    outcube = outds.GetData()
    if propvar:
        invar  = inds.GetVar()
        outvar = outds.GetVar()

    # Instantiate a mapping object between two datasets and pre-calculate
    # the reverse transforms so geometric_transform doesn't have to do
    # everything in a Python loop:
    mapper = astro_ds.PixMapper(inds, outds)
    mapper.all_reverse()

    # Resample the cube:
    ndi.geometric_transform(incube, mapper.invert, outds.shape, outcube,
        order=3, mode='constant', cval=0.0, prefilter=True)

    if propvar:
            ndi.geometric_transform(invar, mapper.invert, outds.shape, outvar,
                order=3, mode='constant', cval=0.0, prefilter=True)

    # Allow correcting flux for the new binning if the data haven't been flux
    # calibrated yet, since ndimage resamples flux density rather than dividing
    # up a fixed number of counts (this is unlikely to be useful often but I
    # did it before realizing that; currently this step assumes the same
    # wavelength binning over all spatial pixels; it also assumes that the
    # wavelength was linear to begin with and if it wasn't will apply the same
    # correction erroneously)
    if fluxcons:

        # Evaluate wavelength at each slice first:
        # (transform method currently doesn't accept co-ords along only the
        # wavelength axis)
        nw = outds.shape[outds.dispaxis]
        idxarr = numpy.zeros((outds.ndim, nw), dtype=numpy.int16)
        idxarr[outds.dispaxis] = numpy.arange(nw)
        warr = outds.TransformToWCS(idxarr)[outds.dispaxis]

        # Get ratio of each increment WRT the reference increment (which is
        # the same as in the equivalent linear case):
        ratio = outds.GetWCSIncrement(waxis=outds.dispaxis, coord=warr) / \
            outds.GetWCSIncrement(waxis=outds.dispaxis)
    
        # Broadcast up the dimensionality of the ratio array to that of the
        # science array and take their product:
        idx = [numpy.newaxis for axis in range(outds.ndim)]
        idx[outds.dispaxis] = slice(None, None, None) # like [:,newaxis,newaxis]
        ratarr = ratio[idx]
        outcube *= ratarr

        # Whether we multiply variance by the same correction factor as the
        # science array or its square depends on whether we're trying to
        # preserve the overall S/N over some wavelength range or the S/N per
        # pixel, respectively (since we're not tracking covariance properly).
        # I think there should be 1*ratarr if the science array is corrected
        # for conservation of flux in e- and another one (whether or not the
        # science array is corrected) if we want to maintain the S/N per pixel,
        # ie. the power of ratarr is somewhere between 0 and 2. Here we choose
        # to preserve overall S/N, as the alternative will be difficult to make
        # sense of later, and do a single correction here (if applicable):
        if propvar: outvar *= ratarr

# End (LogBin function to resample with log wavelength binning)

