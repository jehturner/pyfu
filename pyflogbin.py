# Copyright(c) 2013 Association of Universities for Research in Astronomy, Inc.
# and James E.H. Turner.
#
# 'pyflogbin' Python module for rebinning cubes logarithmically in wavelength
#
# Version  Nov, 2012  JT Initial version

import pyfits, numpy
from scipy import ndimage
import astro_ds

reload(astro_ds)

# Pyflogbin main routine (non-PyRAF interface):
def pyflogbin(inimage, outimage):

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
    outds.SetGridFrom([inds])

    # Do stuff here

    # Copy the SCI header from the input file as a starting point
    # for the output header:
    scihdr = inhdulist['sci',1].header.copy()

    # Append the output DataSet as an extension of the output HDUList:
    outds.WriteHDU(outhdulist, header=scihdr)

    # Close files, writing the output to disk:
    outhdulist.close(output_verify="warn")
    inhdulist.close(output_verify="warn")

# End (pyflogbin main Python routine)


