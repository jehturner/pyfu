# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
# and James Turner
#
# Version  Feb-Apr, 2006  JT Initial test version
#
"""
A module for high-level manipulation of astronomical datasets (especially
IFU data and Gemini-style MEF files) with STScI's numpy and PyFITS
"""
import math, pyfits, numpy, string, pyfu_transform
import numpy.linalg, numpy.numarray.mlab


# Extension naming convention:
sciname = 'SCI'
varname = 'VAR'
dqname  = 'DQ'


# Class describing a single science data array (with associated WCS,
# variance, data quality etc).
class DataSet:
    """Class describing a science dataset, consisting of a data array
    and WCS parameters etc., with optional variance / data quality"""

    # Initialize the DataSet dimensions and WCS to None until they are
    # defined by the appropriate method (eg. _ReadHdr() or SetGridFrom()):
    ndim  = None
    shape = None
    cridx = None
    cd    = None
    icd   = None
    crval = None

    # Initialize file mapping and modified data array to None until they
    # are updated:
    _hdulist  = None
    _filename = None
    _extver   = None
    _newdata  = None
    _newDQ    = None

    # Cop out and implement DQ etc. later on :-). Actually partly done now.

    # Initialize DataSet by checking for specified extension(s) and
    # reading the WCS:
#    def __init__(self, filename, extver=1, create=False, from_ds=None):
    def __init__(self, hdulist=None, extver=1):

        # If we were passed a FITS HDUList, create the DataSet from,
        # and mapped to, that file:
        if isinstance(hdulist, pyfits.HDUList):

            # Keep a record of the HDUList/file mapping
            # (is it safe to use _file attribute?):
            self._hdulist  = hdulist
            self._filename = hdulist[0]._file.name
            self._extver   = extver

            # Make string corresponding to the SCI header name in IRAF:
            scihdrname = self._filename+'['+str(sciname)+','+str(extver)+']'

            # Try to get the specified SCI header from the HDUList:
            try:
                scihdr = hdulist[sciname, extver].header
            except KeyError:
                raise KeyError('extension '+scihdrname+' not found')

            # Parse WCS etc from the SCI extension header:
            self._ReadHdr(scihdr)

        # Error for now if hdulist is not an HDUList (but maybe support
        # filenames later, opening and caching the HDUList):
        elif hdulist is not None:
            raise TypeError('DataSet.__init__(): bad hdulist argument, '\
                            +'%s' % type(hdulist))
            
        # Otherwise, if we're creating a new empty DataSet, don't
        # do any dynamic initialization for now

    # End (DataSet.__init__() method)


    # Private method to read an existing FITS header:
    def _ReadHdr(self, scihdr):

        # It would be good to put in a check for skew in the CD matrix,
        # since all the calculations here assume there isn't any. We
        # should really return an error to the user if things are skewy.

        # Get the number of dimensions (required FITS keyword):
        self.ndim = _ndim = scihdr['naxis']

        # For each image dimension (backwards WRT FITS convention), get the
        # length (required FITS keyword) and reference index (WRT 0):
        dimlist, ridxlist = [], []
        for dim in range(_ndim, 0, -1):

            # Get dimensions:
            dimlist.append(scihdr['naxis'+str(dim)])

            # Get reference index = CRPIX-1:
            keyw = 'CRPIX'+str(dim)
            try:
                cridx = scihdr[keyw]-1.0
            except KeyError:
                # print self._key_warning(keyw, scihdrname, 1)
                cridx = 0.0
            ridxlist.append(cridx)

        # End (loop over input dimensions)

        # Get the CD matrix, assuming for the time being that the number
        # of WCS dimensions equals the number of image dimensions:
        # cdmatrix = [range(_ndim) for i in range(_ndim)] # old: use list
        cdmatrix = numpy.zeros((_ndim,_ndim), numpy.float32)
        for wdim in range (_ndim):
            for dim in range (_ndim):
                keyw = 'CD'+str(_ndim-wdim)+'_'+str(_ndim-dim)
                try:
                    cdval = scihdr[keyw]
                except KeyError:
                    if dim == wdim: cdval = 1.0
                    else: cdval = 0.0
                    # print self._key_warning(keyw, scihdrname, cdval)
                cdmatrix[dim,wdim] = cdval
        # End (loop over CD matrix)

        # For each WCS dimension, get the reference co-ordinate (CRVAL):
        rvallist = []
        for wdim in range(_ndim, 0, -1):
            keyw = 'CRVAL'+str(wdim)
            try:
                crval = scihdr[keyw]
            except KeyError:
                crval = 0.0
                # print self._key_warning(keyw, scihdrname, crval)
            rvallist.append(crval)
        # End (loop over CRVAL keywords)

        # Store the WCS values as Python tuples or lists, as appropriate:
        # (tuples can be passed directly to some functions eg. for 'shape',
        # whereas lists are easily converted to arrays, eg. cd matrix):
        self.shape  = tuple(dimlist)
        self.cridx  = ridxlist
        self.cd     = cdmatrix
        self.crval  = rvallist

        # Should we return a warning if any of the defaults were used
        # because keywords could not be found?
        
    # End (private method to read existing FITS header)


    # Method to map a DataSet in memory to an HDUList extension:
    def WriteHDU(self, hdulist, header=None, extver=1):

    # For the time being, don't keep a record of the mapping to the output
    # HDUList because it may not be written to disk yet and therefore can't
    # delete/read back the data attribute etc. - doesn't behave the same as
    # an HDUList that has been read from file.

#         # Raise error if the DataSet is already mapped to a different file
#         # and do nothing if it is already mapped to the same file:
#         if not self._hdulist==None:
#             if hdulist==self._hdulist: return
#             else:
#                 raise RuntimeError('mapto(): DataSet already mapped to ' \
#                     +str(self._filename))

#         # Store a record of the new mapping:
#         self._hdulist = hdulist
#         if not hdulist[0]._file==None:   # (may not be defined yet)
#             self._filename = hdulist[0]._file.name
#         self._extver = extver

        # When we append a new extension to an HDUList, PyFITS doesn't set
        # the _extver attribute, so we can't subsequently look up HDUs by
        # extension name & version. Thus we can't reliably check whether the
        # extver we want to update/append already exists by its name & number.
        # Loop over the existing HDUs & check the names instead:
        scicount = 0
        scihdu = None
        for hdu in hdulist:
            if isinstance(hdu, pyfits.ImageHDU) and hdu.name==sciname:
                scicount+=1
                try:
                    # In case extver has been set in the header but not
                    # the HDU and there are skipped extvers or what not,
                    # synchronize the count with the header EXTVER:
                    sciver = hdu.header['extver']
                    scicount = sciver
                except KeyError:
                    pass
                if scicount==extver:
                    scihdu = hdu
                    break

        # If the output extver doesn't already exist, append it to the
        # HDUList as a new extension:
        if scihdu==None:

            scihdu = pyfits.ImageHDU(header=header, name=sciname, \
                         data=self.GetData())
            
            hdulist.append(scihdu)

        # If the output extver exists already, just replace the header
        # and data (may not work if old dimensionality doesn't match?):
        else:

            # Use the specified header (may be the same as previously):
            scihdu.header = header

            # scihdu.data = self._newdata  # this only worked by accident
            
            # Update the data array unless the DataSet is already mapped
            # from this same HDU and the data haven't changed:
            if self._newdata is not None \
               or self._hdulist[sciname, self._extver] is not scihdu:
                scihdu.data = self.GetData()

        # Ensure the extension version is correct in the header:
        scihdu.header.update('extver', extver)
                
        # Update SCI header with the WCS etc:
        self._UpdateHdr(scihdu.header)

    # End (method to write a modified dataset)
    

    # Method to update an HDU header with the DataSet WCS etc.
    def _UpdateHdr(self, scihdr):

        # Set dimensionality:
        ndim = self.ndim
        scihdr.update('naxis', ndim)

        ndimplus1 = ndim+1

        # Loop over the FITS dimensions and set the ref. pix/value for each:
        for fdim in range(1, ndimplus1):

            # Get the array dimension, counting backwards WRT FITS:
            adim = ndim-fdim

            # Don't set the axis lengths because that will be done when
            # we write the array data (setting explicitly causes probs)
            #
            # # Set axis length in pixels:
            # keyw = 'NAXIS'+str(fdim)
            # scihdr.update(keyw, self.shape[adim])

            # Set reference pixel for this axis:
            keyw = 'CRPIX'+str(fdim)
            scihdr.update(keyw, self.cridx[adim]+1.0)

            # Set world reference value for this axis:
            keyw = 'CRVAL'+str(fdim)
            scihdr.update(keyw, self.crval[adim])

        # End (loop over FITS dimensions)

        # Loop over the permutations of FITS dimensions and set the
        # CD matrix elements:
        for wdim in range(1, ndimplus1):
            for fdim in range(1, ndimplus1):
                keyw = 'CD'+str(wdim)+'_'+str(fdim)
                scihdr.update(keyw, self.cd[ndim-fdim,ndim-wdim])
        # End (loop over CD matrix)
            
    # End (method to update an HDU header with the DataSet WCS etc)
    

    # Define a warning string for missing header keywords:
    def _key_warning(self, keyword, imname, default):
        warning = 'Warning: no \''+str(keyword)+'\' in '+ \
            str(imname)+': using value '+str(default)
        return warning
    # End (warning string method)


    # Method to return a copy of the DataSet's data array:
    def GetData(self):

        # If the DataSet already has a data array from a previous call,
        # return that:
        if not (self._newdata is None):
            return self._newdata

        # If the DataSet is newly created and the dimensions have been set,
        # return an empty array of the correct size:
        if self._hdulist is None:

            # If the array dimensions haven't yet been specified for the
            # new dataset, return an error (we don't have any information
            # here from which to determine the required size):
            if self.shape==None:
                raise RuntimeError('GetData(): can\'t create a new array' \
                                   + 'before setting its dimensions')

            # Get a new array full of zeros and save a reference to it, so
            # we can save changes to it later:
            self._newdata = darr = numpy.zeros(self.shape, 'float32')
            # Is Float32 the best type to use? Defaults to Int if unspecified!

        # If the DataSet is from an existing file, read the array data:
        else:

            # How can we check whether the HDUList has been closed by the
            # main program? Is there an attribute of HDUList that we can
            # use, or do we have to try operating on the data? The following
            # assignment will work whether the data are readable or not, but
            # operating on the array may subsequently fail in the main
            # program if the file has been closed (eg. if memory mapped).

            # (After figuring this out, we could support re-opening the
            # filename if the HDUList has been closed and then closing the
            # HDUList when the DataSet is deleted or unmapped)

            # Get the data array from the HDUList (then reset it to being
            # lazy in the HDUList so it doesn't persist in memory longer
            # than needed by the calling program):
            try:
                darr = numpy.float32(self._hdulist[sciname, self._extver].data)
                del self._hdulist[sciname, self._extver].data

            # Return an error if the data can't be read because there is
            # no image in the extension:
            except IndexError:
                scihdrname = self._filename+'['+sciname+','+ \
                             str(self._extver)+']'
                raise IndexError('no data in '+scihdrname)

        # Return the read in or newly created data array:
        return darr
        
    # End (method GetData)


    # Method to return a copy of the DataSet's DQ array. This is more of a
    # quick hack for now than for SCI, assuming DQ is read-only.
    def GetDQ(self):

        # If the DataSet already has a DQ array from a previous call,
        # return that:
        if self._newDQ is not None:
            return self._newDQ

        # If the DataSet is newly created and the dimensions have been set,
        # return an empty array of the correct size:
        if self._hdulist is None:

            # If the array dimensions haven't yet been specified for the
            # new dataset, return an error (we don't have any information
            # here from which to determine the required size):
            if self.shape==None:
                raise RuntimeError('GetDQ(): can\'t create a new array' \
                                   + 'before setting its dimensions')

            # Get a new array full of zeros and save a reference to it, so
            # we can save changes to it later:
            self._newDQ = darr = numpy.zeros(self.shape, 'int16')

        # If the DataSet is from an existing file, read the array data:
        else:

            # Get the data array from the HDUList (then reset it to being
            # lazy in the HDUList so it doesn't persist in memory longer
            # than needed by the calling program):
            try:
                darr = numpy.int16(self._hdulist[dqname, self._extver].data)
                del self._hdulist[dqname, self._extver].data

            # Return an error if the data can't be read because there is
            # no image in the extension:
            except IndexError:
                dqhdrname = self._filename+'['+dqname+','+ \
                            str(self._extver)+']'
                raise IndexError('no data in '+dqhdrname)

        # Return the read in array:
        return darr
        
    # End (method GetDQ)


    # Method to get a 2D image by summing the DataSet's data array over
    # any higher dimensions. If the 'match' parameter specifies a reference
    # DataSet, the output image is transformed to its co-ordinate system.
    def GetTelImage(self, match=None):

        # Would be useful to add a check that self and match contain CD
        # matrices of the correct form (rank), to catch programming errors

        # Collapse the array down to a 2D image, assuming the last 2
        # dimensions (ie. the first 2 in IRAF/FITS) are spatial axes:
        image = self.GetData()
        
        for dim in range(0, self.ndim-2):
            image = numpy.numarray.mlab.mean(image,0)

        # If there is a reference dataset, transform image accordingly:
        if match:

            # Get the CD matrix and the inverse of the CD matrix for
            # the reference dataset:
            scd = self.GetCD()
            micd = match.GetICD()

            # Get the 2x2 spatial section of each CD matrix:
            d1,d2 = self.ndim-2, self.ndim
            scd2d = scd[d1:d2,d1:d2]
            d1,d2 = match.ndim-2, match.ndim
            micd2d = micd[d1:d2,d1:d2]

            # Calculate the product of the CD matrix with the inverse
            # of the reference CD matrix, to get the transform.:
            trmat = numpy.dot(scd2d, micd2d)

            # Transform the image onto the appropriately sized array:
            image = pyfu_transform.CentredAffine(image, trmat, cval=0.0)

        return image

    # End (method to get a 2D image from an ND dataset)


    # Method to decide whether DataSet is on the standard celestial
    # orientation (with East anticlockwise from North), based on the WCS.
    def IsStdOrient(self):

        # Assume for now (horror!) that the first 2 FITS axes are spatial
        # and that the first WCS axis is parallel to RA and the second to
        # Dec (also see GetPA function below)

        # Anything 'funny' in the CD matrix (other than a normal rotation,
        # flip and scaling) will produce a False result without further ado.

        # Check the relative signs of the CD matrix coefficients, to find
        # the sense of the WCS axes is relative to the data:
        cdmatrix = self.GetCD()
        fc1, fc2 = self.ndim-1, self.ndim-2
        diag1 = cdmatrix[fc1,fc1] * cdmatrix[fc2,fc2]
        diag2 = cdmatrix[fc1,fc2] * cdmatrix[fc2,fc1]

        if diag1 < 0.0 or diag2 > 0.0: stdorient = True
        else: stdorient = False

        return stdorient
    
    # End (method to test for std orientation)


    # Method to calculate DataSet's PA from its CD matrix
    def GetPA(self):

        # Assume for now that FITS WCS axis 1 is parallel to RA and FITS
        # WCS axis 2 is parallel to Dec (NB. this does NOT assume that the
        # actual *array* axes are parallel to RA and Dec, just the *WCS*)

        # I think this only works for an orthogonal WCS ... translation,
        # rotation and separate scalings in each dimension are allowed,
        # but skew (rotation of one WCS axis with respect to the other)
        # is not. I'm not quite 100% sure, however, that skew isn't the
        # same thing as a rotation, followed by a change of aspect ratio
        # followed by a rotation back, followed by an aspect scaling back.
        
        # Everything works as long as the CD matrix is formed like this
        # (which is what I assumed in the calculations):
        #
        #                              CD1_1             CD1_2
        #       Std orientation: -delt1*cos(theta) -delt2*sin(theta)
        #                        -delt1*sin(theta)  delt2*cos(theta)
        #
        #       N->E clockwise:   delt1*cos(theta) -delt2*sin(theta)
        #                         delt1*sin(theta)  delt2*cos(theta)
        #                              CD2_1             CD2_2
        #
        # ... where theta = -PA is the angle from y to North and
        #     delt1/delt2 are the pixel scales

        # NB. Calculating the angle of the RA and Dec axes for annotation
        # is easier than this, because we can just use CD2_1/CD2_2 and
        # CD1_1/CD1_2. Here, however, we need to factor out delt1 and delt2,
        # to get theta rather than the angle(s) with respect to the grid.

        cdmatrix = self.GetCD()
        fc1, fc2 = self.ndim-1, self.ndim-2
        
        diag1 = cdmatrix[fc1,fc1] * cdmatrix[fc2,fc2]
        diag2 = cdmatrix[fc1,fc2] * cdmatrix[fc2,fc1]

        # After hacking around some trig identities, we get the following
        # value for the magnitude of the position angle for (-90<theta<=90)
        # (this is definitely correct):
        try:
            mag_theta = 90.*math.acos((diag1 + diag2) / (diag1 - diag2)) \
                         /math.pi
        except (ZeroDivisionError, ValueError):
            # Something is 'wrong' with the CD matrix, eg. the WCS axes
            # are parallel, so the PA calculation fails
            return None

        # If WCS axes 1 & 2 are reversed, we would insead get the following
        # (with a sign flip in the argument of arccos):
        # mag_theta = 90.*math.acos((diag1 + diag2) / (diag2 - diag1))/math.pi

        # The above 2 results do not depend on whether there is a sign flip
        # in of one of the array axes WRT N up E left

        # Now we need to adjust the quadrant according to whether
        # abs(theta)<90 or 90<abs(theta)<180; use the dDec/dy term to check,
        # since it is always positive up to 90 and then negative:
        if cdmatrix[fc2,fc2] < 0.0: mag_theta = 180. - mag_theta

        # Get the sign of the PA from the dRA/dy term, which always has
        # the same sign as the PA:
        if cdmatrix[fc1,fc2] < 0.0: theta = -mag_theta
        else: theta = mag_theta

        return theta

    # End (method to calculate PA from CD matrix)


    # Method to calculate DataSet's scale factors from its CD matrix
    def GetScale(self):

        # Here we assume that the first 2 FITS axes are spatial and
        # that there is no rotation that would mix spatial and spectral
        # co-ordinates in a single axis. We also (I think) assume no skew,
        # as for GetPA. Thus we are assuming the following form for the
        # CD matrix (it doesn't make any difference here if East is
        # clockwise from North instead of on the standard orientation):
        #
        #   delt3                 0                 0
        #       0   delt2*cos(theta) -delt1*sin(theta)
        #       0  -delt2*sin(theta) -delt1*cos(theta)
 
        cdmatrix = self.GetCD()
        fc1, fc2 = self.ndim-1, self.ndim-2

        deltlist = []

        # For dimensions higher than 2, just append the diagonal CD
        # elements to the list of scales:
        for dim in range(self.ndim-2):
            deltlist.append(cdmatrix[dim,dim])

        # Calculate the spatial scale(s), allowing for rotation between
        # the 2 spatial axes. The reported spatial scales are always
        # positive, so to reproduce the CD matrix fully, GetScale() and
        # GetPA() need to be used in conjunction with IsStdOrient(), as
        # follows:
        #
        # if IsStdOrient(): sgn = -1 else: sgn = 1
        # 
        # CD2_2 = delt2*math.cos(math.radians(-pa))
        # CD2_1 = sgn*delt1*math.sin(math.radians(-pa))
        # CD1_2 = -delt2*math.sin(math.radians(-pa))
        # CD1_1 = sgn*delt1*math.cos(math.radians(-pa))

        delt = math.sqrt(cdmatrix[fc2,fc1]**2 + cdmatrix[fc2,fc2]**2)
        deltlist.append(delt)

        delt = math.sqrt(cdmatrix[fc1,fc1]**2 + cdmatrix[fc1,fc2]**2)
        deltlist.append(delt)

        return tuple(deltlist)

    # End (method to calculate scale factors from the CD matrix)


    # Method to set the WCS parameters and array size based on a list of
    # input datasets, optionally overriding specific WCS parameters :
    def SetGridFrom(self, dslist, pa=None, stdorient=None, scale=None):

        # Like other DataSet methods, this currently assumes (only when
        # overriding WCS parameters) that the first 2 axes are spatial and
        # there is no rotation between spatial and spectral dimensions
        # Later use a dictionary to identify which are the spatial axes?

        # Set output dimensionality based on the first input dataset:
        ds1 = dslist[0]
        ndim = self.ndim = ds1.ndim

        # Set output CD matrix based on the first input dataset,
        # overriding any WCS parameters specified as arguments:
        if pa==None and stdorient==None and scale==None:

            # If no WCS parameters are overriden, just copy the input
            # CD matrix and its inverse, if present:
            self.cd = ds1.GetCD().copy()
            self.icd = numpy.array(ds1.GetICD(force=False))
            if self.icd == numpy.array(None): self.icd=None
                
        else:
            
            # If at least one WCS parameter is overriden, recalculate
            # the CD matrix after deriving the other parameters from the
            # input CD matrix.

            # Add more type checking later, to catch silly program errors?

            # Get PA from input if not defined (as an int or float)
            if pa==None: pa = ds1.GetPA()

            # Get flip from input if not defined (as True/False)
            if stdorient==None: stdorient = ds1.IsStdOrient()
            if stdorient: sign = -1
            else: sign = 1

            # Get scale(s) from input if not defined (as tuple of floats)
            if scale==None:  
                scale = ds1.GetScale()
            elif not type(scale)==tuple:
                raise TypeError('SetGridFrom(): \'scale\' must be a tuple')
            elif not len(scale)==ndim:
                raise ValueError('SetGridFrom(): \'scale\' should match ' +\
                                 'number of WCS dimensions')

            # (Later on, could allow spatial scales only to be specified,
            # treating a length-2 tuple as applying to the spatial axes
            # and a single float as applying to both spatial axes?)

            # Start with a zeroed CD matrix of the right dimensionality:
            # cdmatrix = [[0.0 for j in range(ndim)] for i in range(ndim)]
            cdmatrix = numpy.zeros((ndim,ndim), numpy.float32)

            # Figure out which are the spatial dimensions (assume for
            # now that they are the first 2 FITS axes):
            d1,d2 = ndim-2, ndim-1

            # Loop over the non-spatial diagonal terms of the CD matrix
            # and enter the appropriate scale values:
            for dim in range(d1):
                cdmatrix[dim,dim] = scale[dim]

            # Calculate the spatial terms of the CD matrix from the scale,
            # PA & sign (remember d1==FITS_axis_2 & d2==FITS_axis_1):
            pa_rad = math.radians(pa)
            cdmatrix[d2,d2] =  sign * scale[d2] * math.cos(pa_rad)
            cdmatrix[d1,d2] =         scale[d1] * math.sin(pa_rad)
            cdmatrix[d2,d1] = -sign * scale[d2] * math.sin(pa_rad)
            cdmatrix[d1,d1] =         scale[d1] * math.cos(pa_rad)

            # Store the final CD matrix, setting the inverse matrix to
            # None for now, since TransformFromWCS will calculate it below:
            self.cd = cdmatrix
            self.icd = None

        # End (if WCS parameters are not overriden)

        # Loop over the input datasets and get the corners bounding each
        # dataset, relative to the output system:
        allcorners = []
        for dataset in dslist:

            # Calculate the dataset's array corners:
            corners = pyfu_transform.GetCorners(dataset.shape)

            # Transform the corners to world co-ords, common to all datasets:
            corners = [dataset.TransformToWCS(corner) for corner in corners]

            # Transform WCS corners to relative co-ords on the output grid
            # (with an arbitrary zero point):
            corners = [self.TransformFromWCS(corner, relative=True) \
                       for corner in corners]

            # Append this datset's corners to the overall list:
            allcorners.extend(corners)

        # End (loop over datasets, calculating corners)

        # Calculate the required axis lengths from the min & max corner
        # co-ordinate along each axis, keeping a copy of the minima:
        outshape, mindex = [], []
        for axis in range(ndim):
            cvals = [corner[axis] for corner in allcorners]
            mincval = min(cvals)
            mindex.append(mincval)
            outshape.append(int(max(cvals)-mincval+1))

        # Store the required output array size:
        self.shape = tuple(outshape)

        # For the time being, set the reference index to zero along
        # each axis (might prefer to set it to the middle of the cube
        # spatially, but avoids assuming which axes are spatial here)
        self.cridx = [0 for axis in range(ndim)]

        # Calculate the WCS zero points at index (0, 0, ... 0) by
        # transforming the minimum pixel indices in the above relative
        # output pixel co-ordinates back to WCS co-ordinates:
        self.crval = list(self.TransformToWCS(mindex, relative=True))

    # End (method to set the WCS parameters & array size)


    # Method to transform array co-ordinates or offsets to world
    # co-ordinates:
    def TransformToWCS(self, coords, relative=False):

        # If the co-ordinate tuple is shorter than the CD matrix, match
        # the final rows/columns (ie. the first FITS axes):
        d1,d2 = self.ndim-len(coords),self.ndim

        # For absolute coords, subtract reference coordinate (CRPIX-1):
        if not relative:
            cridx = self.cridx[d1:d2]
            coords = [coord-ref for coord,ref in zip(coords,cridx)]

        # Multiply the tuple by the CD matrix:
        cd = self.GetCD()[d1:d2,d1:d2]
        coords = numpy.dot(cd, coords)

        # For absolute co-ords, add the WCS zero point (CRVAL):
        if not relative:
            crval = self.crval[d1:d2]
            coords = [coord+ref for coord,ref in zip(coords,crval)]

        return tuple(coords)

    # End (method to transform array co-ordinates to world co-ordinates)


    # Method to transform array co-ordinates or offsets from world
    # co-ordinates to array (pixel-1) co-ordinates:
    def TransformFromWCS(self, coords, relative=False):

        # If the co-ordinate tuple is shorter than the CD matrix, match
        # the final rows/columns (ie. the first FITS axes):
        d1,d2 = self.ndim-len(coords),self.ndim

        # For absolute co-ords, subtract the WCS zero point (CRVAL):
        if not relative:
            crval = self.crval[d1:d2]
            coords = [coord-ref for coord,ref in zip(coords,crval)]

        # Multiply the tuple by the inverse CD matrix:
        icd = numpy.array(self.GetICD())[d1:d2,d1:d2]
        coords = numpy.dot(icd, coords)

        # For absolute coords, add the reference coordinate (CRPIX-1):
        if not relative:
            cridx = self.cridx[d1:d2]
            coords = [coord+ref for coord,ref in zip(coords,cridx)]

        return tuple(coords)
        
    # End (method to transform world co-ordinates to array co-ordinates)
    

    # Method to apply an adjustment to a DataSet's WCS:
    def OffsetWCS(self, offset):

        # Error if the offset argument isn't a tuple or list:
        if not (type(offset)==tuple or type(offset)==list):
            raise TypeError('OffsetWCS(): offset must be a tuple or list')

        # If the offset list/tuple is shorter than the number of
        # dimensions, match with the first FITS axes:
        d1,d2 = self.ndim-len(offset), self.ndim

        # Check that we didn't get too many dimensions in the offset:
        if d1 < 0:
            raise ValueError('OffsetWCS(): more offset components than WCS'
                             +' dimensions')

        # Add the offsets to the crval offsets:
        crvala = self.crval[0:d1]
        crvalb = [rval+oval for rval,oval in zip(self.crval[d1:d2], offset)]

        # Store the new reference co-ordinates (=CRVALs):
        self.crval = crvala + crvalb

    # End (method to apply an offset to a DataSet's WCS)


    # Method to return a DataSet's CD matrix as an array:
    def GetCD(self):

        if not isinstance(self.cd, numpy.ndarray):
            raise ValueError('GetCD(): CD matrix undefined for \'' \
                             + self._filename+'\'')

        return self.cd
    
    # End (method to return a DataSet's CD matrix)


    # Method to return a DataSet's inverse CD matrix as an array:
    def GetICD(self, force=True):

        if self.icd is None and force:
            if not isinstance(self.cd, numpy.ndarray): #self.cd==None:
                raise ValueError('GetICD(): CD matrix undefined for \'' \
                                 + self._filename+'\'')
            else:
                self.icd = numpy.linalg.inv(self.cd)

        return self.icd
    
    # End (method to return a DataSet's inverse CD matrix)


# End (class DataSet)


# Function to open a list of filenames as a list of HDUList objects:
def OpenMEFList(inlist, mode='readonly'):
    """Open a list of files as a list of PyFITS HDUList objects"""

    # To support list files and wildcards, used irafglob to expand out the
    # file list before passing to this function.

    # Loop over the input files:
    meflist = []
    for filename in inlist:

        # Attempt to open the current file:
        try:
            hdulist = pyfits.open(filename, mode=mode)

        # If it fails, close any previously opened files before re-posting
        # the exception:
        except IOError:
            for hdulist in meflist:
                hdulist.close()
            raise

        # If it worked, add the file to the list:
        meflist.append(hdulist)

    # Return the list of open PyFITS HDULists:
    return meflist

# End (function to open file list as HDUList objects)


# Function to close a list of PyFITS HDUList objects:
def CloseMEFList(inlist):
    """Close a list of PyFITS HDUList objects"""

    # Loop over the input files:
    for hdulist in inlist:
        hdulist.close(output_verify="warn")

# End (function to close a list of HDUList objects)


# Function to ensure a filename includes a '.fits' extension:
def WithFITS(filename):

    if not (filename.endswith('.fits') or filename.endswith('.fit')):
        filename = filename + '.fits'

    return filename

# End (function to include '.fits' extension)


# Function to 'extract' a list of DataSets from one or more HDULists:
def GetDataSets(inlist, extver=1):

    """Open list of DataSets from a list of PyFITS HDULists or filenames"""

    # If input is a single MEF file, convert to a list:
    if isinstance(inlist, pyfits.HDUList) or isinstance(inlist, str):
        inlist = [inlist]

    # If the extension version is a single integer, convert to a list:
    if isinstance(extver, int):
        extver = [extver]

    # If the extension version is a string other than '*', return an
    # error for the time being (could parse ranges in future):
    if isinstance(extver, str):
        if not extver=='*':
            raise TypeError('GetDataSets(): bad extver argument \'%s\'' \
                            % extver)

    # Loop over the input HDULists (or filenames):
    dslist = []
    for hdulist in inlist:

        # If extver is a wildcard, get a list of SCI extensions:
        if extver=='*': extvlist = SCIExtVer(hdulist)
        else: extvlist = extver

        # Loop over the SCI extensions:
        for ver in extvlist:

            # Old behaviour opened the HDU and passed the filename
            # instead of the HDUList

            # Create a new DataSet from the current extension
            # (pass the HDUList and specify which HDU):
            ds = DataSet(hdulist, ver)

            # Append dataset to the list:
            dslist.append(ds)

        # End (loop over the SCI extensions)

    return dslist
      
# End (function to extract a list of DataSets from HDULists)


# Function to return the SCI ext version numbers for an HDUList:
def SCIExtVer(hdulist):

    verlist = []

    # Currently no checking is done for duplicate extension versions
    # (including blank extver, which implies 1)
    
    # Check that the input is an hdulist:
    if not isinstance(hdulist, pyfits.HDUList):
        raise TypeError('SCIExtVer(): argument must be an HDUList, not ' + \
                        '\'%s\'' % type(hdulist))

    # Loop over the SCI extensions:
    for hdu in hdulist:
        if string.upper(hdu.name)==sciname:

            # Get the version number for the current ext:
            try:
                ver = int(hdu.header['EXTVER'])
            except:
                ver = 1

            # Append the version number to the list:
            verlist.append(ver)

    # Later on add support for simple FITS files, returning extver=0?

    # Return the list of extension versions:
    return verlist

# End (function to return the SCI ext versions for an HDUList)


# # Function to convert a list of image names to a list of pyifu DataSet
# # objects (without opening the data arrays, to avoid keeping everything in
# # memory before needed):
# def DataSetList(inlist):
#     """Create a list of astro_ds.DataSet objects from an image list"""

#     # Loop over the input files:
#     dslist = []
#     for imname in inlist:

#         # Get the dimensions, WCS etc:
#         try:
#             datasetobj = DataSet(imname)
#         except IOError:
#             # Don't do anything special for just now
#             raise
#         dslist.append(datasetobj)

#     return dslist

# # End (function to create a list of DataSet objects)

