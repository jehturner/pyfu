# Copyright(c) 2006 Association of Universities for Research in Astronomy, Inc.
#
# Version  Feb-Apr, 2006  JT Initial test version
#
#   Python library modules for the pyfu Pyraf package, to help with
#   N-dimensional image transformations,

import numpy
import numpy.linalg
from scipy import ndimage

# Function to call numpy's affine transformation with the output
# automatically recentred & resized:
def CentredAffine(input, matrix, offset=0.0, output_type=None, order=3,
                  mode='constant', cval=0.0, prefilter=True):

    # Get the size & dimensionality of the input array:
    inshape = numpy.array(input).shape
    naxis = len(inshape)

    # Invert the transformation matrix:
    invmatrix = numpy.linalg.inv(matrix)

    # Generate a list of corner co-ordinates for the input array:
    corners = GetCorners(inshape)

    # Transform the corners to the output system:
    corners = [tuple(numpy.dot(invmatrix, corner)) for corner in corners]

    # Calculate the length of each output axis from the min & max co-ords
    # of the corners:
    outshape = []
    for axis in range(naxis):
        cvals = [corner[axis] for corner in corners]
        outshape.append(int(max(cvals)-min(cvals)+1))
    outshape = tuple(outshape)

    # Calculate the offset needed to keep the transformed data centred in
    # the array:
    incen  = [0.5*(axlen-1) for axlen in inshape]
    outcen = [0.5*(axlen-1) for axlen in outshape]
    cenoff = numpy.array(incen) - numpy.dot(matrix, outcen)

    # Transform the data onto the new array:
    output = ndimage.affine_transform(input, matrix,
             offset=offset+cenoff, output_shape=outshape,
             output=output_type, order=order, mode=mode, cval=cval,
             prefilter=prefilter)

    return output

# End (shifted & resized affine transformation function)

    
# Recursive function to calculate the corner indices of an array of the
# specified shape:
def GetCorners(shape):

    if not type(shape)==tuple:
        raise TypeError('GetCorners argument is non-tuple')

    if len(shape)==1:
        corners = [(0,), (shape[0]-1,)]
    else:
        shape_less1    = shape[1:len(shape)]
        corners_less1  = GetCorners(shape_less1)
        corners = []
        for corner in corners_less1:
            newcorner = (0,) + corner
            corners.append(newcorner)
            newcorner = (shape[0]-1,) + corner
            corners.append(newcorner)
        
    return corners

# End (function to calculate the corners of an array)

