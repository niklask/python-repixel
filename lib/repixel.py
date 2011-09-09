#
# Python/numarray repixeling package
#  module: repixel.py
#

import sys
import pickle
from numarray import *

from mathx import *
import geometry

# --------------------------------------------------------------------------------
#  Function definitions
# --------------------------------------------------------------------------------

# Function repixelCartesianToCartesian2D
#
# Perform a pixelized transformation of a 2D array
#
# Changelog:
# 08/11/05: Function renamed
# 08/11/05: Changed the way coordinates are given by the user. Now the center of
#           the grids/arrays is given by the coordinates of its center. The rotation
#           is done in the local coordinate system of grid/array B and then
#           positioned in the global coordinate system.
#
# Basic idea:
#  For each pixel of array B do the following
#  1. find the corner vertices
#  2. find the pixels of A that are partially or fully overlapped
#  3. for each such pixel, calculate area of intersection
#  4. value of pixel is the sum of all parts in 3 normalized to the pixel area
#
def repixelCartesianToCartesian2D(A, B, posA, posB, theta,
                                  pxlsizeA=array([1.0,1.0]), pxlsizeB=array([1.0,1.0])):

    # calculate the coordinates of all pixel corners
    # matrix Bx is such that Bx[i,j] returns the x coordinate of the upper vertex of pixel (i,j)
    # of array B; By, Ax, Ay works the same way

    index = indices((B.shape[0]+1, B.shape[1]+1), type=Float32)

    # first position the grid relative to the origin of the _local_ system
    cenB = B.shape*pxlsizeB/2.
    Bx = index[0,:,:]*pxlsizeB[0] - cenB[0]
    By = index[1,:,:]*pxlsizeB[1] - cenB[1]
    # then perform a rotation in 2D about the _local_ origin
    if (theta != 0.0):
        Bxp, Byp = geometry.rotate2D([Bx, By], rad(theta))
    else:
        Bxp = Bx
        Byp = By
    # finally position the grid/array in _global_ system using posB
    # (posB gives the coordinates for the center (i.e. origin)
    Bxp = Bxp + posB[0]
    Byp = Byp + posB[1]

    index = indices((A.shape[0]+1, A.shape[1]+1), type=Float32)
    cenA = A.shape*pxlsizeA/2

    Ax = index[0,:,:]*pxlsizeA[0] - cenA[0] + posA[0]
    Ay = index[1,:,:]*pxlsizeA[0] - cenA[1] + posA[1]

    # store the boundaries of the A array
    Ax_min = Ax.min()
    Ax_max = Ax.max()
    Ay_min = Ay.min()
    Ay_max = Ay.max()

    # time to calculate the values of pixels in array B
    # this is done pixel by pixel in two nested for-loops
    n = B.shape
    i_range = range(n[0])
    j_range = range(n[1])

    i_max = A.shape[0] - 1
    j_max = A.shape[1] - 1

    for j in j_range:
        for i in i_range:
            # find the four corner vertices of pixel (i,j)
            # each column represent a point and row 0 is x and row 1 is y
            pixel = array([[Bxp[i+1, j], Bxp[i+1, j+1], Bxp[i, j+1], Bxp[i, j]],
                           [Byp[i+1, j], Byp[i+1, j+1], Byp[i, j+1], Byp[i, j]]])

            # find the lowest and highest x,y for the points above
            x_min = pixel[0,:].min()
            x_max = pixel[0,:].max()
            y_min = pixel[1,:].min()
            y_max = pixel[1,:].max()

            # if the pixel is entirely outside the A array then assign a value 0
            if ((x_max < Ax_min) or (x_min > Ax_max) or
                (y_max < Ay_min) or (y_min > Ay_max)):
                B[i,j] = 0.0
            # otherwise we must calculate the average
            else:
                # find the indices of overlapped pixels of array A
                k_min = int(floor((x_min + cenA[0])/pxlsizeA[0]))
                k_max = int(floor((x_max + cenA[0])/pxlsizeA[0]))
                l_min = int(floor((y_min + cenA[1])/pxlsizeA[1]))
                l_max = int(floor((y_max + cenA[1])/pxlsizeA[1]))

                k_range = k_min + arange((k_max - k_min) + 1)
                l_range = l_min + arange((l_max - l_min) + 1)

                # this is the pixel of B we want to calculate value for
                B_pixel = [pixel[:,0], pixel[:,1], pixel[:,2], pixel[:,3]]

                # init value to 0.0
                value = 0.0
                # now loop over all the overlapped pixels
                for l in l_range:
                    # do not process if this is a pixel outside the boundaries
                    if ((l < 0) or (l > j_max)):
                        continue
                    
                    for k in k_range:
                        # do not process if this is a pixel outside the boundaries
                        if ((k < 0) or (k > i_max)):
                            continue

                        # if the underlying pixel has value close to zero then
                        # there is no point in calculating intersection
                        if (fcmp(A[k,l], 0.0)):
                            continue

                        # this is the pixel of A we want to calculate intersection with 
                        p1 = array([Ax[k+1, l], Ay[k+1, l]])
                        p2 = array([Ax[k+1, l+1], Ay[k+1, l+1]])
                        p3 = array([Ax[k, l+1], Ay[k, l+1]])
                        p4 = array([Ax[k, l], Ay[k, l]])
                        # represent it as a polygon, i.e. a list of vertices
                        A_pixel = [p1, p2, p3, p4]

                        # calculate the area of intersection between A_pixel and B_pixel
                        intersection = geometry.intersectPolyPoly2D(A_pixel, B_pixel)
                        if (intersection is not None):
                            A_area = geometry.polygonArea2D(A_pixel)
                            i_area = geometry.polygonArea2D(intersection)
                            value += (A[k,l] * i_area/A_area)

                B[i,j] = value

    return

#
# Function repixel3DCartesianToCartesian3D
#
# Perform a pixelized transformation of a 3D array
#
# Changelog:
# 08/11/05: Function renamed
#
# Basic idea:
#  For each pixel of array B do the following
#  1. find the corner vertices => calculate normal vectors, 6 sides, and the 12 side vectors
#  2. find the pixels of A that are partially or fully overlapped
#  3. for each such pixel, calculate volume of intersection
#  4. value of pixel is the sum of all parts in 3 normalized to the pixel volume
#  (i.e. essentially the same as the 2D version except that we have to deal with polyhedrons
#   and volumes instead of polygons and areas)
#
def repixelCartesianToCartesian3D(A, B, posA, posB, phi, theta, psi,
                                  pxlsizeA=array([1.0, 1.0, 1.0]), pxlsizeB=array([1.0, 1.0, 1.0])):

    # calculate the coordinates of all pixel corners
    # matrix Bx is such that Bx[i,j,k] returns the x coordinate of the upper vertex of pixel (i,j,k)
    # of array B; By, Bz, Ax, Ay, Az works the same way

    index = indices((B.shape[0]+1, B.shape[1]+1, B.shape[2]+1), type=Float32)

    # first position the grid relative to the origin of the _local_ system
    cenB = B.shape*pxlsizeB/2.
    Bx = index[0,:,:,:]*pxlsizeB[0] - cenB[0]
    By = index[1,:,:,:]*pxlsizeB[1] - cenB[1]
    Bz = index[2,:,:,:]*pxlsizeB[2] - cenB[2]

    # 3D rotation around the _local_ origin (using Euler angles)
    Bxp, Byp, Bzp = geometry.rotate3D([Bx, By, Bz], rad(phi), rad(theta), rad(psi))

    # finally position the grid/array in _global_ system using posB
    # (posB gives the coordinates for the center (i.e. origin)
    Bxp = Bxp + posB[0]
    Byp = Byp + posB[1]
    Bzp = Bzp + posB[2]
  
    index = indices((A.shape[0]+1, A.shape[1]+1, A.shape[2]+1), type=Float32)
    cenA = A.shape*pxlsizeA/2.
    Ax = index[0,:,:,:]*pxlsizeA[0] - cenA[0] + posA[0]
    Ay = index[1,:,:,:]*pxlsizeA[1] - cenA[1] + posA[1]
    Az = index[2,:,:,:]*pxlsizeA[2] - cenA[2] + posA[2]

    # store the boundaries of the A array
    Ax_min = Ax.min()
    Ax_max = Ax.max()
    Ay_min = Ay.min()
    Ay_max = Ay.max()
    Az_min = Az.min()
    Az_max = Az.max()

    # create the normal vectors of each array
    # there are six normal vectors, one for each face and every pixel has
    # the same set of normal vectors
    N = [array([1.0, 0.0, 0.0]), array([0.0, 1.0, 0.0]), array([0.0, 0.0, -1.0]),
         array([-1.0, 0.0, 0.0]), array([0.0, -1.0, 0.0]), array([0.0, 0.0, 1.0])]
    Np = []
    for i in range(len(N)):
        Npx, Npy, Npz = geometry.rotate3D(N[i], rad(phi), rad(theta), rad(phi))
        Np.append([Npx, Npy, Npz])

    # time to calculate the values of pixels in array B
    # this is done pixel by pixel in three nested for-loops
    n = B.shape
    k_range = range(n[2])
    j_range = range(n[1])
    i_range = range(n[0])

    i_max = A.shape[0] - 1
    j_max = A.shape[1] - 1
    k_max = A.shape[2] - 1

    for k in k_range:
        for j in j_range:
            for i in i_range:
                # find the vertices of the pixel (i,j,k), there are 8
                # each column is a vertex and the rows are x,y and z
                pixel = array([[Bxp[i+1,j+1,k], Bxp[i+1,j+1,k+1], Bxp[i,j+1,k+1], Bxp[i,j+1,k], Bxp[i+1,j,k], Bxp[i+1,j,k+1], Bxp[i,j,k+1], Bxp[i,j,k]],
                               [Byp[i+1,j+1,k], Byp[i+1,j+1,k+1], Byp[i,j+1,k+1], Byp[i,j+1,k], Byp[i+1,j,k], Byp[i+1,j,k+1], Byp[i,j,k+1], Byp[i,j,k]],
                               [Bzp[i+1,j+1,k], Bzp[i+1,j+1,k+1], Bzp[i,j+1,k+1], Bzp[i,j+1,k], Bzp[i+1,j,k], Bzp[i+1,j,k+1], Bzp[i,j,k+1], Bzp[i,j,k]]])

                # find the lowest and highest x,y,z for the points above
                x_min = pixel[0,:].min()
                x_max = pixel[0,:].max()
                y_min = pixel[1,:].min()
                y_max = pixel[1,:].max()
                z_min = pixel[2,:].min()
                z_max = pixel[2,:].max()

                # if the pixel is entirely outside the A array then assign a value 0
                if ((x_max < Ax_min) or (x_min > Ax_max) or
                    (y_max < Ay_min) or (y_min > Ay_max) or
                    (z_max < Az_min) or (z_min > Az_max)):
                    B[i,j,k] = 0.0
                # otherwise we must calculate the average
                else:
                    # find the indices of overlapped pixels of array A
                    l_min = int(floor((x_min + cenA[0])/pxlsizeA[0]))
                    l_max = int(floor((x_max + cenA[0])/pxlsizeA[0]))
                    m_min = int(floor((y_min + cenA[1])/pxlsizeA[1]))
                    m_max = int(floor((y_max + cenA[1])/pxlsizeA[1]))
                    n_min = int(floor((z_min + cenA[2])/pxlsizeA[2]))
                    n_max = int(floor((z_max + cenA[2])/pxlsizeA[2]))

                    l_range = l_min + arange((l_max - l_min) + 1)
                    m_range = m_min + arange((m_max - m_min) + 1)
                    n_range = n_min + arange((n_max - n_min) + 1)

                    # this is the pixel of B we want to calculate the value for
                    # the sides are 12 segments binding together the verticies
                    B_pixel_S = [[pixel[:,0], pixel[:,1]], [pixel[:,1], pixel[:,2]],
                                 [pixel[:,2], pixel[:,3]], [pixel[:,3], pixel[:,0]],
                                 [pixel[:,4], pixel[:,5]], [pixel[:,5], pixel[:,6]],
                                 [pixel[:,6], pixel[:,7]], [pixel[:,7], pixel[:,4]],
                                 [pixel[:,0], pixel[:,4]], [pixel[:,1], pixel[:,5]],
                                 [pixel[:,2], pixel[:,6]], [pixel[:,3], pixel[:,7]]]
                    # six facets
                    B_pixel_V = [pixel[:,0], pixel[:,0], pixel[:,0],
                                 pixel[:,6], pixel[:,6], pixel[:,6]]
                        
                    # init value to 0.0
                    value = 0.0
                    # now loop over all overlapped pixels
                    for n in n_range:
                        # do not process if this is a pixel outside the boundaries
                        if ((n < 0) or (n > k_max)):
                            continue

                        for m in m_range:
                            # do not process if this is a pixel outside the boundaries
                            if ((m < 0) or (m > j_max)):
                                continue

                            for l in l_range:
                                # do not process if this is a pixel outside the boundaries
                                if ((l < 0) or (l > i_max)):
                                    continue

                                # if the underlying pixel has value close to zero then
                                # there is no point in calculating intersection
                                if (fcmp(A[l,m,n], 0.0)):
                                    continue

                                # this is the pixel of A we want to calculate intersection with
                                p1 = array([Ax[l+1,m+1,n], Ay[l+1,m+1,n], Az[l+1,m+1,n]])
                                p2 = array([Ax[l+1,m+1,n+1], Ay[l+1,m+1,n+1], Az[l+1,m+1,n+1]])
                                p3 = array([Ax[l,m+1,n+1], Ay[l,m+1,n+1], Az[l,m+1,n+1]])
                                p4 = array([Ax[l,m+1,n], Ay[l,m+1,n], Az[l,m+1,n]])
                                p5 = array([Ax[l+1,m,n], Ay[l+1,m,n], Az[l+1,m,n]])
                                p6 = array([Ax[l+1,m,n+1], Ay[l+1,m,n+1], Az[l+1,m,n+1]])
                                p7 = array([Ax[l,m,n+1], Ay[l,m,n+1], Az[l,m,n+1]])
                                p8 = array([Ax[l,m,n], Ay[l,m,n], Az[l,m,n]])

                                # represent it as a set of segments and facets
                                A_pixel_S = [[p1, p2], [p2, p3], [p3, p4], [p4, p1],
                                             [p5, p6], [p6, p7], [p7, p8], [p8, p5],
                                             [p1, p5], [p2, p6], [p3, p7], [p4, p8]]
                                A_pixel_V = [p1, p1, p1, p7, p7, p7]

                                # calculate the volume of intersection
                                intersection = geometry.intersectPolyhedra3D(A_pixel_S, A_pixel_V, N,
                                                                             B_pixel_S, B_pixel_V, Np)
                                if (intersection is not None):
                                    i_volume = geometry.polyhedronVolume3D(intersection)
                                    A_volume = geometry.polyhedronVolume3D(A_pixel_V)
                                    if (i_volume is not None):
                                        value += (A[l,m,n] * i_volume/A_volume)
                                    
                    B[i,j,k] = value

    return

#
# Function repixelSphereToSphere
#
# Perform a pixelized transformation of a 2D array defined on the surface of a sphere
# That is, considering coordinates (l,b) instead of (x,y)
#
# Changelog:
#  08/31/05: Made it possible to have different pixel sizes on the A and B grids
#            B array is now sent as argument, just like the other repixeling functions
#
# ToDo:
#  * calculate values at polar regions
#
# Basic idea:
#  Transform each pixel and calculate its value as the sum of the pixels of the underlying
#  coordinate system, normalized to the area
#
#  Coordinates of A are (l,b) and coordinates of B are (l',b'). For the repixelation, (l',b')
#  are mapped to (l,b) with an inverse rotation (three Euler angles). This rotation is such
#  that e_{x'} lies in the direction of point A, given in (l,b), and the north pole is on the great
#  circle (still in (l,b)).
#
# TODO: fix repixeling at the transition l=-180 -> l=180
#
def repixelSphereToSphere(A, B, pointB, pointA=array([0.0, 0.0])):
    # 11/17/06: removed pixel size defs from function definition
    #           since coordinate system never changes we can calculate pixel sizes
    #           from the array shapes
    pxlsizeA = array([180.0/A.shape[0], 360.0/A.shape[1]])
    pxlsizeB = array([180.0/B.shape[0], 360.0/B.shape[1]])
    
    # calculate the coordinates of all pixel corners
    # matrix Bglon is such that Bglon[i,j] returns the glon coordinate of the upper vertex of pixel (i,j)
    # of array B; Bglat, Aglon, Aglat works the same way
    index = indices((A.shape[0]+1, A.shape[1]+1), type=Float32)

    # the coordinate system for A is defined such that glon=0, glat=0 points toward the
    # galactic center, thus, rows of the array will represent constant glat, starting with
    # glat=90 and decreasing, and columns represent constant glon, extending from -180 to 180
    Aglat = 90.0 - index[0,:,:]*pxlsizeA[1]
    Aglon = index[1,:,:]*pxlsizeA[0] - 180.0
    Ax,Ay,Az = geometry.convertGCToCartesian(Aglon, Aglat)

    # and we calculate Bglon,Bglat,Bx,By,Bz in the same manner as above
    Bglat = 90.0 - index[0,:,:]*pxlsizeB[1]
    Bglon = index[1,:,:]*pxlsizeB[0] - 180.0
    Bx,By,Bz = geometry.convertGCToCartesian(Bglon, Bglat)

    # then rotate these coordinates
    # phi, theta, psi are the Euler angles for the rigid body rotation as defined in
    # http://mathworld.wolfram.com/EulerAngles.html
    # with this definition, unit vector e_{x} is rotated to e_{x'} which points through
    # coordinates l,b given by pointB

    dl = pointB[0] - pointA[0]
    db = pointB[1] - pointA[1]
    # phi: rotate CW about z-axis, here e_{x'}->e_{y} 
    phi = -90.0
    # theta: rotate CW about x-axis, tilt e_{x'} up/down
    # ISSUE: in EGRET data b=-90->90 and not 90->-90 as defined above
    #        this means that theta here has opposite sign
    theta = db
    # psi: rotate CW about z-axis again, line up e_{x'} toward pointB
    psi = 90.0 - dl
    
    # rotate all the grid points in x,y,z space using the above Euler angles
    Bx,By,Bz = geometry.rotate3D([Bx,By,Bz], rad(phi), rad(theta), rad(psi))
    # and convert back to galactic coordinates l,b
    # (this will overwrite the Bglon and Bglat from above, but we do not need those)
    Bglon,Bglat = geometry.convertCartesianToGC(Bx, By, Bz)

    # time to calculate the values of the pixels in the transformed array B
    # i will loop over constant glon in B
    # j will loop over constant glat in B
    # for now, we are excluding the polar regions, i.e. 90 -> 90-dglat and -90 -> -90+dglat
    n = B.shape
    i_range = range(1, n[0] - 1)
    j_range = range(n[1])

    cenA = [A.shape[0]*pxlsizeA[1]/2., A.shape[1]*pxlsizeA[0]/2.]

    i_max = A.shape[0]

    # loop over all all pixels (elements) of array B and calculate the values
    # the first and last row is excluded because they represent the polar cap, treat them
    # separately in the end
    # j loops over rows and i loops over columns in the array; 0,0 is the first elements
    for j in j_range:
        for i in i_range:
            # find the four vertices that make up the pixel (i,j)
            # each column represent a point and row 0 is glon and row 1 is glat
            pixel = array([[Bglon[i+1,j], Bglon[i+1,j+1], Bglon[i, j+1], Bglon[i, j]],
                           [Bglat[i+1,j], Bglat[i+1,j+1], Bglat[i, j+1], Bglat[i, j]]])
                
            # find the boundary glon and glat for this pixel
            glon_min = pixel[0,:].min()
            glon_max = pixel[0,:].max()
            glat_min = pixel[1,:].min()
            glat_max = pixel[1,:].max()

            k_min = int(floor((cenA[0] - glat_max)/pxlsizeA[1]))
            k_max = int(floor((cenA[0] - glat_min)/pxlsizeA[1]))
            l_min = int(floor((cenA[1] + glon_min)/pxlsizeA[0]))
            l_max = int(floor((cenA[1] + glon_max)/pxlsizeA[0]))

            k_range = k_min + arange((k_max - k_min) + 1)
            l_range = (l_min + arange((l_max - l_min) + 1))%A.shape[1]

            B_pixel = [pixel[:,0], pixel[:,1], pixel[:,2], pixel[:,3]]

            l_max = A.shape[1]
                
            # init B pixel value to 0.0
            value = 0.0
            
            # now loop over all the overlapped pixels in A
            for l in l_range:
                # the break flag is set to 1 when the pixel is found to overlap with the
                # polar region
                breakflag = 0
                
                for k in k_range:
                    # do not process this pixel if it is part of the polar region, i.e.
                    # the first and last rows of pixels in A
                    if ((k < 0) or (k >= i_max)):
                        breakflag = 1
                        break

                    # if the underlying pixel has value close to zero then
                    # there is no point in calculating intersection
                    if (fcmp(A[k,l], 0.0)):
                        continue
                                
                    # this is the pixel of A we want to calculate intersection with
                    p1 = array([Aglon[k+1, l], Aglat[k+1, l]])
                    p2 = array([Aglon[k+1, (l+1)%l_max], Aglat[k+1, (l+1)%l_max]])
                    p3 = array([Aglon[k, (l+1)%l_max], Aglat[k, (l+1)%l_max]])
                    p4 = array([Aglon[k, l], Aglat[k, l]])
                    # represent it as a polygon, i.e. a list of vertices
                    A_pixel = [p1, p2, p3, p4]

                    # calculate the area of intersection between A_pixel and B_pixel
                    intersection = geometry.intersectPolyPoly2D(A_pixel, B_pixel)
                    if (intersection is not None):
                        i_area = geometry.polygonArea2D(intersection)
                        A_area = geometry.polygonArea2D(A_pixel)
                        if (i_area < 0.0) or (A_area < 0.0):
                            continue
                        value += (A[k,l] * i_area/A_area)
                        
                if (breakflag == 1):
                    # if we are breaking because of overlap with polar region then reset
                    # pixel value to 0.0
                    # TODO: when the polar region is overlapped, we should do better than this
                    value = 0.0
                    break

            B[i,j] = value
            
#
# Test functions
#
# The functions below are used to test the functionality of the repixeling functions
#
def test1_repixel2D():
    print "Test 1: repixeling a 2D grid"

    # define the array we want to repixelize
    # the ones make up a triangle
    A = array([[1.0, 1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0],
               [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0]])
    # and the one we want to repixelize into
    B = zeros((10,10), type=Float32)

    # the coordinate of the grids/arrays relative to the origin in the
    # _global_ system
    posA = array([0.0, 0.0])
    posB = array([0.0, 0.0])

    # print the original array
    print A
    
    # do the repixelization
    # the coordinates defined and an angle theta is a complete description
    # of the transformation (rotation + translation)
    repixelCartesianToCartesian2D(A, B, posA, posB, 90.0, pxlsizeB=array([0.5, 0.5]))

    # print the result
    print B

def test2_repixel3D():
    print "Test 2: repixeling a 3D grid"

    # define the array we want to repixelize
    # the ones make up a pyramid
    A = array([[[1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
               [[1.0, 1.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 0.0]],
               [[1.0, 1.0, 1.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]])

    # and the one we want to repixelize into
    B = zeros((3,3,3), type=Float32)

    # the coordinate of the upper left corner of each array
    posA = array([0.0, 0.0, 0.0])
    posB = array([0.0, 0.0, 0.0])


    print A
    # do the repixelization
    # the coordinates defined and three Euler angles is a complete description
    # of the transformation (rotation + translation)
    # phi = -90, theta=90, psi=90 rotates the grid so that the axes goes: x' = x, z' = y, y' = -z
    # (above angles might not give mentioned rotation since we had the Euler angles incorrectly defined)
    repixelCartesianToCartesian3D(A, B, posA, posB, 0.0, 90.0, 0.0)

    print B

def test3_repixelSphereToSphere():
    print "Test 3: repixeling of a 2D grid on a sphere"
    print "        input and output arrays are pickled to files"
    print "        input: sphere_A.num; output: sphere_B.num"
    
    # define an array to map on the sphere
    A = zeros((36, 72), type=Float32)

    # create a 6 pixel wide line around the equator
    m = A.shape[0]
    n = A.shape[1]
    A[m/2-3:m/2+3,0:n] = ones((6,n), type=Float32)

    # and the one we want to repixelize into
    B = zeros((36, 72), type=Float32)

    f = open("sphere_A.num", "w")
    pickle.dump(A, f)
    f.close()

    # define the point through which the new galactic center should point
    pointB = array([90.0, 0.0])

    # do the repixelization
    # the coordinates defined gives the Euler angles for the rotation of the sphere
    repixelSphereToSphere(A, B, pointB)

    f = open("sphere_B.num", "w")
    pickle.dump(B, f)
    f.close()

#
# Function runTest
#
# Main entry point for the test functions. Takes what test to run as argument.
#
def runTest(testno = 0):
    if (testno == 1):
        # Test 1: repixel of 2D grid
        test1_repixel2D()
    elif (testno == 2):
        # Test 2: repixel of 3D grid
        test2_repixel3D()
    elif (testno == 3):
        # Test 3: repixel of a 2D grid on a sphere
        test3_repixelSphereToSphere()

#
# Function main
#
def main(argv = None):
    if (argv is None):
        argv = sys.argv
    argc = len(argv)

    if (argc > 1):
        test = int(argv[1])
        runTest(test)

if __name__ == "__main__":
    sys.exit(main())
    
