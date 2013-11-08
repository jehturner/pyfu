# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
#
# PyRAF interface functions for pyfu.pyfmosaic and related routines
#
# Version  Feb-May, 2006  JT Initial test version

from pyraf import iraf
import pyfmosaic, irafglob_jt, astro_ds

# Temporary - reload module development changes:
reload(pyfmosaic)
reload(irafglob_jt)


# Pyfmosaic PyRAF interface:
#def pyfmosaic_iraf(inimages, outimage, posangle, method, cenfirst, cenlast):
def pyfmosaic_iraf(inimages, outimage, posangle, separate):

    # Convert inimages string to a Python list of filenames using a
    # modified version of STScI irafglob that returns an error when any
    # specified files are missing and expands out ".fits" extensions
    inlist = irafglob_jt.irafglob(inimages, matchfits=True, onfail='error')
    if len(inlist)==0:
        raise IOError('no files matching \'%s\'' % inimages)

    # Ensure the output filename includes a '.fits' extension:
    outimage = astro_ds.WithFITS(outimage)

    # Get output PA (if INDEF, use first input file's WCS):
    if posangle == iraf.INDEF: posangle = None

    # Call the main Python routine:
    pyfmosaic.pyfmosaic(inlist, outimage, posangle, separate)

# End (pyfmosaic PyRAF interface routine)


# Pyfalign PyRAF interface:
def pyfalign_iraf(images, method):

    # Convert inimages string to a Python list of filenames using a
    # modified version of STScI irafglob that returns an error when any
    # specified files are missing and expands out ".fits" extensions
    inlist = irafglob_jt.irafglob(images, matchfits=True, onfail='error')
    if len(inlist)==0:
        raise IOError('no files matching \'%s\'' % images)

    # Call the main Python routine:
    pyfmosaic.pyfalign(inlist, method)

# End (pyfalign PyRAF interface routine)


# Define pyfmosaic as an IRAF task:
pyfmosaic_parfile = iraf.osfn("pyfu$pyfmosaic.par")
tsk = iraf.IrafTaskFactory(taskname="pyfmosaic", value=pyfmosaic_parfile,
        function=pyfmosaic_iraf)

# Define pyfalign as an IRAF task:
pyfalign_parfile = iraf.osfn("pyfu$pyfalign.par")
tsk = iraf.IrafTaskFactory(taskname="pyfalign", value=pyfalign_parfile,
        function=pyfalign_iraf)

