# Copyright(c) 2013 Association of Universities for Research in Astronomy, Inc.
#
# PyRAF interface functions for pyfu.pyflogbin and related routines
#
# Version  Nov, 2013  JT Initial version
# Version  Nov, 2018  JT Python 3 compatibility

from pyraf import iraf

from pyfu import pyflogbin, astro_ds


# Pyflogbin PyRAF interface:
def pyflogbin_iraf(inimage, outimage, flux, var):

    # Ensure the input filename includes a '.fits' extension:
    inimage = astro_ds.WithFITS(inimage)

    # Ensure the output filename includes a '.fits' extension:
    outimage = astro_ds.WithFITS(outimage)

    # Call the main Python routine:
    pyflogbin(inimage, outimage, fluxcons=flux, propvar=var)

# End (pyflogbin PyRAF interface routine)


pyflogbin_parfile = iraf.osfn("pyfu$pyflogbin.par")
tsk = iraf.IrafTaskFactory(taskname="pyflogbin", value=pyflogbin_parfile,
        function=pyflogbin_iraf)

