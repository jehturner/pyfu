# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
# and James Turner
#
# 'pyfmosaic' main Python module for mosaicing IFU datacubes
#
# Version  Feb-May, 2006  JT Initial test version
#              May, 2011  JT Version with basic support for input DQ
#              Oct, 2013  JT Fix ndimage import & finally in version control!

import pyfits, numpy
from scipy import ndimage
from stsci import imagestats
from astropy.nddata.convolution import convolve as acnv
import astro_ds
#import numdisplay

# Temporary - reload module development changes:
reload(astro_ds)


# Pyfmosaic main routine (non-PyRAF interface):
#
# list   inimages = Python list of input image names
#                   (no wildcards here, only in the PyRAF interface)
# str    outimage = output image name
# float  posangle = output position angle
#
def pyfmosaic(inimages, outimage, posangle=None):

    """Mosaic IFU datacubes, based on WCS"""

    # Add a bit more sanity checking for the input parameters when
    # the task is closer to being finished?

    # Open a list of input files as PyFITS HDUList objects:
    meflist = astro_ds.OpenMEFList(inimages)

    # Get a list of DataSets from the input files:
    dslist = astro_ds.GetDataSets(meflist)

    # Copy the output PHU from the first input file:
    outphu = pyfits.PrimaryHDU(header=meflist[0][0].header.copy())
    # outhdulist = pyfits.HDUList(outphu) # old: now re-open instead

    # Save the PHU to disk, to ensure the file is writeable:
    outphu.writeto(outimage)
    del outphu

    # Re-open the minimal output image for appending:
    outhdulist = pyfits.open(outimage, mode='append')

    # Do the alignment in separate pyfalign task instead of here:
    # # Measure a list of offsets from the list of input datasets and
    # # update their WCS offset accordingly:
    # if not method=='wcs': AdjOffsets(dslist)

    # Create an empty DataSet for the output and set the WCS and array
    # size, based on the WCS and array sizes of the list of input datasets:
    outds = astro_ds.DataSet()
    outds.SetGridFrom(dslist, pa=posangle)

    # Interpolate input cubes onto the output:
    AddCubes(dslist, outds)

    # Copy the SCI header from the first input file as a starting point
    # for the output header:
    scihdr = meflist[0]['sci',1].header.copy()

    # Append the output DataSet as an extension of the output HDUList:
    # (should we maintain a copy of the full SCI header in each DataSet
    # instead of passing it here?)
    outds.WriteHDU(outhdulist, header=scihdr)

    # Write the output file to disk & close it:
    outhdulist.close(output_verify="warn")

    # Close the input files:
    astro_ds.CloseMEFList(meflist)

# End (pyfmosaic main Python routine)


# Pyfalign main routine (non-PyRAF interface):
#
# list     images = Python list of image names to update with new WCS
#                   (no wildcards here, only in the PyRAF interface)
# str      method = datacube alignment method
#
def pyfalign(images, method):

    """Spatially align WCS of IFU datacubes onto a common system"""

    known_methods=['centroid', 'correlate']

    # Only one alignment method is implemented for just now:
    if method not in known_methods:
        raise 'pyfalign: invalid alignment method: \'%s\'' % method

    # Open a list of input files as PyFITS HDUList objects:
    meflist = astro_ds.OpenMEFList(images, mode='update')

    # Get a list of DataSets from the input files:
    dslist = astro_ds.GetDataSets(meflist)

    # Measure a list of offsets from the list of input datasets and
    # update their WCS offset accordingly:
    AdjOffsets(dslist, method)

    # Update each input HDUList with the corresponding modified DataSet:
    for hdulist, ds in zip(meflist, dslist):
        ds.WriteHDU(hdulist, header=hdulist['sci',1].header)

    # Close the input files:
    astro_ds.CloseMEFList(meflist)

# End (pyfalign main Python routine)


# Function to derive spatial offsets for a list of datacube arrays:
def AdjOffsets(dslist, method='centroid'):
    """Measure spatial offsets for a list of >=2D datasets and adjust
       their WCS parameters to put them all on the same co-ordinates"""

    # Loop over the input datasets:
    for dataset in dslist:

        # Get an image summed over wavelength and nominally transformed to
        # the co-ordinate system of the first DataSet:
        image = dataset.GetTelImage(match=dslist[0])

        if method=='centroid':

            # Measure the centroid of the brightest source:
            pkcen = ImagePeakCen(image)

            # (use a more sophisticated algorithm later, but don't assume that
            # all peaks stay within the field at all pointings)

            # Transform the peak co-ordinates to RA/Dec:
            pkwcs = dataset.TransformToWCS(pkcen)
            if dataset==dslist[0]: pkref = pkwcs

        elif method=='correlate':

            # Unlike centroid, this method currently requires that the input
            # cubes have the same PA & scale, as no transformation is
            # performed on the pixel data before cross-correlating. It should
            # work with images of different sizes (untested).

            # Remember the reference image and the World co-ordinates of its
            # centre, to which we'll apply measured cross-correlation offsets
            # to get the corrected centre co-ordinates of the other images:
            if dataset==dslist[0]:
                refim = image
                # pkoff = dataset.TransformToWCS([0 for axis in image.shape])
                pkwcs = dataset.TransformToWCS([0.5*(axlen-1) for axlen \
                    in image.shape])
                pkref = pkwcs

            # For subsequent datasets, find the sub-pixel shift from upsampled
            # cross-correlation, then call TransformToWCS to convert to World
            # offsets: 
            else:
                pixoff = ImageCorrelationShifts(refim, image)
                pkwcs = dataset.TransformToWCS([0.5*(axlen-1)-axoff for \
                    axlen, axoff in zip(image.shape, pixoff)])

        # Calculate the adjustment needed to match the first dataset's WCS:
        pkoff = [ref-pk for ref,pk in zip(pkref, pkwcs)]

        # Update this dataset's WCS to remove the calculated difference:
        dataset.OffsetWCS(pkoff)

# End (function to derive spatial offsets from list of cube arrays)


# Function for measuring the brightest peak centroid in an image:
def ImagePeakCen(image):
    """Find the centroid of the brightest peak in an image"""

    # This is just a fairly quick-and-dirty algorithm for now

    # Detect all sources at >2*sigma in the image:
    bgstats = imagestats.ImageStats(image,nclip=3,lsig=3.0,usig=2.0)
    imgmask = image >= bgstats.mean + 2*bgstats.stddev
    src_labels, src_num = ndimage.label(imgmask)

    # Identify the label corresponding to the brightest source:
    peaks = numpy.array(ndimage.maximum(image, src_labels, range(1, src_num+1)))
    maxsource = numpy.argmax(peaks)+1   # peaks has maxima for labels 1 onwards

    # Centroid on the area of the source > half-max:
    peaksize = peaks[maxsource-1]-bgstats.mean
    src_labels = (image >= bgstats.mean + 0.5*peaksize) * src_labels
    objcen = ndimage.center_of_mass(image, numpy.float32(src_labels),
               maxsource)  # Why does src_labels need converting to float?
                           # (needed when porting to Numpy)

    return objcen

# End (function for measuring brightest peak centroid)


# Function for finding image shifts via cross-correlation:
def ImageCorrelationShifts(image1, image2):
    """Find the sub-pixel shift between images by cross-correlation"""

    precision = 0.1
    
    # Make a copy of each image with the mean subtracted, to avoid any
    # artificial cross-correlation peak at 0,0.
    mean1 = imagestats.ImageStats(image1, nclip=0).mean
    mean2 = imagestats.ImageStats(image2, nclip=0).mean
    image1 = image1 - mean1
    image2 = image2 - mean2

    # Upsample each image by x10 along each axis, to get 0.1 pixel precision
    # from the discrete cross-correlation without fitting the peak:
    sfactor = 1. / precision
    image1 = ndimage.zoom(image1, sfactor, order=3)
    image2 = ndimage.zoom(image2, sfactor, order=3)

    # Cross-correlate by convolving with the complex conjugate of the
    # reversed reference image:
    # image2 = image1   # DEBUG
    product = acnv.convolve_fft(numpy.conjugate(image1[::-1,::-1]), image2)

    # if dataset==dslist[3]:
    #   pyfits.writeto("im1.fits", image1)
    #   pyfits.writeto("im2.fits", image2)
    #   pyfits.writeto("prod.fits", product)

    # For now we just take the maximum value in the cross-correlation
    # product as the applicable peak. In a test on a single line emission
    # region with asymmetric structure this works well, with a peak that's
    # broad but distinguishable to the nearest upsampled pixel:
    pkoff = numpy.unravel_index(numpy.argmax(product), product.shape)

    # Centre co-ordinates of the output (corresponding to 0 shift):
    pkref = [0.5*(axlen-1) for axlen in product.shape]

    # Shifts WRT the centre:
    pkoff = [precision*(ref-pk) for ref,pk in zip(pkref, pkoff)]

    return pkoff

# End (function for finding image shifts via cross-correlation)


# Function to co-add datacubes with spatial offsets:
def AddCubes(dslist, outds):
    """Co-add DataSets with spatial offsets"""

    # Get an empty output array:
    outcube = outds.GetData()
    # outDQ = outds.GetDQ()
    outshape = outcube.shape

    # Get inverse CD matrix for the output cube:
    outicd = outds.GetICD()

    # Create a mask counting the number of input pixels contributing
    # to each output pixel (maximum 256, no checking!):
    outmask = numpy.zeros(outshape, dtype='uint8')

    # Get world co-ords of the output cube origin:
    outzero = outds.TransformToWCS([0.0 for axis in range(outds.ndim)])

    # Loop over the input datasets:
#    for dataset in dslist[0:1]: # (for testing cube4.fits interpolation)
    for dataset in dslist:

        # Get array data from DataSet:
        cube = dataset.GetData()
        gooddq = 1-dataset.GetDQ()  # Added 2010

        # print 'ref count is ', sys.getrefcount(cube)

        # Calculate the transformation matrix from relative output co-ords
        # to relative input co-ords:
        incd = dataset.GetCD()
        trmatrix = numpy.dot(incd, outicd)

        # Calculate the offset needed to match the WCS zero points:
        troffset = dataset.TransformFromWCS(outzero)

        # Transform the datacube to the output grid:
        trcube = ndimage.affine_transform(cube, trmatrix, troffset,
                 order=3, mode='constant', cval=numpy.nan,
                 output_shape=outshape)

        # Added 2010 to use input DQ now that gfcube is correcting
        # atmospheric dispersion and propagating DQ for the blank edges.
        trdq = ndimage.affine_transform(gooddq, trmatrix, troffset,
               order=1, mode='constant', cval=numpy.nan,
               output_shape=outshape)

        # # Comment old masking, 2010:
        # # Flag blank output pixels, not overlapping this dataset:
        # blank = numpy.isnan(trcube)
        # #
        # # Increment the output mask count at non-blank pixels:
        # outmask += numpy.where(blank, 0, 1)
        
        # # Add the transformed input cube to the output array
        # # (converting any blank pixels to zero):
        # outcube += numpy.where(blank, 0., trcube)

        # New masking using transformed input DQ (haven't thought hard
        # about optimizing the weighting WRT the science array):
        outmask += numpy.where(trdq < 0.5, 0., trdq)
        outcube += numpy.where(trdq < 0.5, 0., trcube)

        # numdisplay.display(outmask[1000,:,:], frame=1)
        # numdisplay.display(outcube[1000,:,:], frame=2)

    # End (loop over input datasets)

    # Save the output mask:
    #outDQ += outmask

    # Replace zeros in the output mask by ones before dividing:
    outmask = numpy.where(outmask, outmask, 1.)

    # Divide co-added output cube by its mask, to convert the sum
    # of pixels to a mean (otherwise different parts of the cube
    # would be normalized differently):
    outcube /= outmask

    # TEST:
#    numdisplay.display(outds.moddata[1000,:,:], frame=1)
#    numdisplay.display(outmask[1000,:,:], frame=2)

# End (function to co-add cubes with spatial offsets)

