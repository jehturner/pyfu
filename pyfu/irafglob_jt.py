"""

License: http://www.stsci.edu/resources/software_hardware/pyraf/LICENSE

"""
import glob

__author__ = 'Paul Barrett; James Turner'
__version__ = '1.0JT'  # (modified from 1.0 by James Turner, May 2006)

def irafglob(inlist, matchfits=False, onfail=None, atfile=None):
    """ Returns a list of filenames based on the type of IRAF input.
 
    Handles lists, wild-card characters, and at-files.  For special
    at-files, use the atfile keyword to process them.

    This function is recursive, so IRAF lists can also contain at-files
    and wild-card characters, e.g. 'a.fits, @file.lst, *flt.fits'.
    """

    # Determine which form of input was provided:
    if isinstance(inlist, list):
        #  python list
        flist = []
        for f in inlist:
            flist += irafglob(f, matchfits=matchfits, onfail=onfail)
    elif ',' in inlist:
        #  comma-separated string list
        flist = []
        for f in inlist.split(','):
            f = f.strip()
            flist += irafglob(f, matchfits=matchfits, onfail=onfail)
    elif inlist[0] == '@':
        #  file list
        flist = []
        for f in open(inlist[1:], 'r').readlines():
            f = f.rstrip()
            # hook for application specific atfiles.
            if atfile:
                f = atfile(f)
            flist += irafglob(f, matchfits=matchfits, onfail=onfail)
    else:
        #  shell globbing
        #  should we remove any image extension/section suffix first?
        flist = glob.glob(inlist)
        #  if there is no exact match, try glob with a '.fits' extension
        #  (should we also match '.fit' and/or other image extensions?)
        if matchfits and len(flist)==0:
            flist = glob.glob(inlist+'.fits')
        #  if there is still no match, optionally raise an exception
        #  unless the input name contains wildcards (in which case a match
        #  is not required; this comparison assumes there are no escaped
        #  wildcards and ignores character ranges, which look like
        #  extensions/sections)
        if len(flist)==0:
            if onfail=='error':
                if not '*' in inlist and not '?' in inlist:
                    raise IOError('cannot open \'%s\'' % inlist)
            elif onfail=='warning':
                if not '*' in inlist and not '?' in inlist:
                    print('Warning: cannot open \'%s\'' % inlist)
            elif onfail=='allerror':
                raise IOError('cannot open \'%s\'' % inlist)
            elif onfail=='allwarning':
                print('Warning: cannot open \'%s\'' % inlist)
            
    return flist

