# Copyright(c) 2013 Association of Universities for Research in Astronomy, Inc.
#
# PyRAF interface functions for pyfu.pyflogbin and related routines
#
# Version  Nov, 2013  JT Initial version

from pyraf import iraf
import pyflogbin, astro_ds

# Temporary - reload module development changes:
reload(pyflogbin)


# Pyflogbin PyRAF interface:
def pyflogbin_iraf(inimage, outimage):

    # Ensure the input filename includes a '.fits' extension:
    inimage = astro_ds.WithFITS(inimage)

    # Ensure the output filename includes a '.fits' extension:
    outimage = astro_ds.WithFITS(outimage)

    # Call the main Python routine:
    pyflogbin.pyflogbin(inimage, outimage)

# End (pyflogbin PyRAF interface routine)


pyflogbin_parfile = iraf.osfn("pyfu$pyflogbin.par")
tsk = iraf.IrafTaskFactory(taskname="pyflogbin", value=pyflogbin_parfile,
        function=pyflogbin_iraf)

