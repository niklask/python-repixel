#
# Python/numarray repixeling package
#   module: geometry.py
#

#
# Basic definitions for geometric entities used within this package:
# 1. nd vertex (point)   numarray with shape (n,)
# 2. vector              basically the same as a vertex
# 3. line segment        list of two vertices representing the end points
# 4. polygon             a list of n vertices ordered CCW, such that
#                          S_i=[V_i, V_i+1], i = 0 -> (n-1), V_n = V_0
#                        is an edge of the polygon
# 5. polyhedron          a list of n normal vectors and a list of n points, such that
#                          {N_i, V_i}
#                        defines a face of the polygon with no specific order between faces
#
# We have chosen not to implement classes for the above mentioned entities. The reason
# for this is to avoid as much overhead over numarray as possible and to maximize the
# use of numarray.
#

from numarray import *

from mathx import *
import hull

#
# Function convertGCToCartesian
#
# l,b given in degrees
#
def convertGCToCartesian(l, b):
    z = sin(rad(b))
    x = sqrt(1 - z**2)*cos(rad(l))
    y = sqrt(1 - z**2)*sin(rad(l))

    return x, y, z
#
# Function convertCartesianToGC
#
# l,b given in degrees
#
def convertCartesianToGC(x, y, z):
    b = arcsin(z)
    l = arctan2(y,x)

    return deg(l), deg(b)

#
# Function rotate2D
#
# Rotate a vector u = (ux, uy) using a rotation matrix
# (can also be applied to pairs of equally sized matrices)
#
def rotate2D(u, theta):
    up_x = u[0]*cos(theta) - u[1]*sin(theta)
    up_y = u[0]*sin(theta) + u[1]*cos(theta)

    return up_x, up_y

#
# Function rotate3D
#
# Rotate a vector u = (ux, uy, uz) using a rotation matrix
# (can also be applied to pairs of equally sized matrices)
# Angles are the Euler angles defined as http://mathworld.wolfram.com/EulerAngles.html
#
def rotate3D(u, phi, theta, psi):
    a11 = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
    a12 = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
    a13 = sin(psi)*sin(theta)
    a21 = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
    a22 = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
    a23 = cos(psi)*sin(theta)
    a31 = sin(theta)*sin(phi)
    a32 = -sin(theta)*cos(phi)
    a33 = cos(theta)

    up_x = a11*u[0] + a12*u[1] + a13*u[2]
    up_y = a21*u[0] + a22*u[1] + a23*u[2]
    up_z = a31*u[0] + a32*u[1] + a33*u[2]

    return up_x, up_y, up_z
#
# Function rotate3Dinverse
#
# Rotate a vector u = (ux, uy, uz) using a rotation matrix
# (can also be applied to pairs of equally sized matrices)
# Angles are the Euler angles defined as http://mathworld.wolfram.com/EulerAngles.html
#
def rotate3Dinverse(up, phi, theta, psi):
    a11 = cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi)
    a12 = -sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi)
    a13 = sin(phi)*sin(theta)
    a21 = cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi)
    a22 = -sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi)
    a23 = -cos(phi)*sin(theta)
    a31 = sin(theta)*sin(psi)
    a32 = sin(theta)*cos(psi)
    a33 = cos(theta)

    u_x = a11*up[0] + a12*up[1] + a13*up[2]
    u_y = a21*up[0] + a22*up[1] + a23*up[2]
    u_z = a31*up[0] + a32*up[1] + a33*up[2]

    return u_x, u_y, u_z

#
# Function planeDistPoint
#
# Distance of a point to a plane
# plane is given by a vertex and a normal vector {V, N}
#
def planeDistPoint(V, N, P):
    sn = -dot(N, (P-V))
    sd = dot(N, N)
    sb = sn / sd
    
    return fabs(sb)

#
# Function triangle2DAarea
#
# Calculate area of a 2D triangle
# (a triangle is a special case of a polygon with only three vertices)
# (taken from: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm)
#
def triangleArea2D(V):
    area = ((V[1][0] - V[0][0])*(V[2][1] - V[0][1]) -
            (V[2][0] - V[0][0])*(V[1][1] - V[0][1])) / 2.0

    return area

#
# Function triangleArea3D
#
# Calculate area of planar 3D triangle
# (a triangle is a special case of a polygon with only three vertices)
#
# NOT IMPLEMENTED!!
#
def triangleArea3D(V):
    return 0.0

#
# Function polygonArea2D
#
# Calculate area of a 2D polygon
# (taken from: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm)
#
def polygonArea2D(V):
    # number of vertices in polygon
    n = len(V)

    area = 0.0
    for i in range(n):
        area += V[(i+1)%n][0] * (V[(i+2)%n][1] - V[i][1])

    return area / 2.0

#
# Function polygonArea3D
#
# Calculate area of planar 3D polygon
# (taken from: http://softsurfe.com/Archive/algorithm_0101/algorithm_0101.htm)
#
def polygonArea3D(V, N):
    # number of vertices in polygon
    n = len(V)

    # select largest abs coordinate to ignore
    a = abs(n)
    coord = 3
    if (a[0] > a[1]):
        if (a[0] > a[2]):
            coord = 1
    elif (a[1] > a[2]):
        coord = 2

    area = 0.0

    # compute area of 2D projection
    for i in range(n):
        if (coord == 1):
            area += V[(i+1)%n][1] * (V[(i+2)%n][2] - V[i][2])
        if (coord == 2):
            area += V[(i+1)%n][0] * (V[(i+2)%n][2] - V[i][2])
        if (coord == 3):
            area += V[(i+1)%n][0] * (V[(i+2)%n][1] - V[i][1])

    # scale to get area before projection
    an = sqrt(a[0]**2 + a[1]**2 + a[2]**2)
    if (coord == 1):
        area *= (an / (2*a[0]))
    if (coord == 2):
        area *= (an / (2*a[1]))
    if (coord == 3):
        area *= (an / (2*a[2]))
    
    return area

#
# Function polygonFromSeg2D
#
# Connect n segments into a list of vertices v representing the polygon
# segments Vi and Vj share a point such that Vi[1] = Vj[0]
#
def polygonFromSeg2D(S):
    # number of segments
    n = len(S)
    p = []

    p.append(S[0][0])
    a = S[0][1]

    for i in range(n-1):
        for j in range(n):
            b = S[j][0]
            if (fcmp(a[0], b[0]) and fcmp(a[1], b[1])):
                a = S[j][1]
                p.append(S[j][0])
                break

        # if the next point is same as the first one added then break
        # this is a crude way of ensuring that two points don't occur twice
        if (fcmp(a[0], p[0][0]) and fcmp(a[1], p[0][1])):
            break

    return p

#
# Function polyhedronCompact3D
#
# Compactifies a set of segments to a set of points, i.e. the vertices of
# the polyhedron
#
def polyhedronCompact3D(S):
    # number of segments
    n = len(S)
    p = []

    p.append(S[0][0])
    p.append(S[0][1])

    for i in range(n-1):
        for j in range(2):
            m = len(p)
            inlist = 0
            for k in range(m):
                a = S[i+1][j]
                b = p[k]

                if (fcmp(a[0], b[0]) and fcmp(a[1], b[1]) and fcmp(a[2], b[2])):
                    inlist = 1
                    break

            if (inlist == 0):
                p.append(S[i+1][j])

    return p

#
# Function tetrahedronVolume3D
#
# Calculate the volume of a tetrahedron, closed polyhedron with four triangular faces
#
def tetrahedronVolume3D(V):
    # number of vertices
    n = len(V)
    if (n != 4):
        return None

    # Note: We need the normal vector to calculate the height and the normal vector
    # is calculated using the cross product. The norm of the normal vector |n|
    # is also directly proportional to the base area. Thus we calculate the area this
    # way to save some computations. Also, we need not explictly calculate the norm.

    u = V[1] - V[0]
    v = V[2] - V[0]
    # n = u x v (cross product)
    n = array([u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]])

    d = fabs(dot(n, V[3] - V[0]))
    
    volume = (1./6.)*d

    return volume

#
# Function polyhedronVolume3D
#
# Calculate the volume of a polyhedron given only by a set of vertices (which are known to form
# a closed polyhedron.)
# Note: if the faces are known then the volume calculation is much easier.
#
def polyhedronVolume3D(V):
    # number of vertices
    n = len(V)

    if (n < 4):
        return None
    elif (n == 4):
        return tetrahedronVolume3D(V)
    else:
        volume = 0.
        
        # find the center of the polyhedron
        Cx = 0.
        Cy = 0.
        Cz = 0.
        for i in range(n):
            Cx += V[i][0]
            Cy += V[i][1]
            Cz += V[i][2]
        C = array([Cx/n, Cy/n, Cz/n])
        
        # calculate the 3D complex hull from the vertices
        # the result is a list of triangles, i.e. each triangle is a list of three
        # numarray vertices
        H = hull.convexHull3D(V)
        n_H = len(H)
        for i in range(n_H):
            t_vol = tetrahedronVolume3D([H[i][0], H[i][1], H[i][2], C])
            volume += t_vol
            
        return volume

    return None

#
# Function cn_PnPoly
#
# Determine inclusion of a point P within a polygon V (n vertices)
# (taken from: http://softsurfer.com/Archive/algorithm_0103/algorithm_0103.htm)
# returns 1 if inside, 0 if not
#
def cn_PnPoly(P, V):
    cn = 0

    n = len(V)
    for i in range(n):
        if (((V[i][1] <= P[1]) and (V[(i+1)%n][1] > P[1])) or
            ((V[i][1] > P[1]) and (V[(i+1)%n][1] <= P[1]))):
            vt = (P[1] - V[i][1]) / (V[(i+1)%n][1] - V[i][1])
            if (P[0] < V[i][0] + vt * (V[(i+1)%n][0] - V[i][0])):
                cn += 1
                
    return (cn & 1)

#
# Function cn_PnPolyhedron
#
# Determine inclusion of a point P within a polyhedron F (n faces)
# the polyhedron is defined by n faces F (the normal vector)
# a point is inside the polyhedron if it lies behind all bounding faces,
# i.e. dot(P-V_i, F_i) <= 0
# as soon as the point lies outside one face, we can stop
# Returns 1 if inside, 0 if not
#
def cn_PnPolyhedron(P, V, F):
    n = len(F)
    for i in range(n):
        d = dot(P-V[i], F[i])
        if (d > 0):
            return 0
        
    return 1

#
# Function intersectSegPoly2D
#
# Find the intersections between a line segment S (two points) and a convex 2D polygon V
# (n vertices)
# (taken from: http://softsurfer.com/Archive/algorithm_0111/algorithm_0111.htm)
# Returns either the segment that intersects or None if error
#
# TODO: should this be modified to work with polygons in 3D as well?
#       cn_PnPoly works in 2D only??
#
def intersectSegPoly2D(S, V):
    # given a 2D vector u, return vector n such that dot(u,n)=0, i.e. orthogonal
    # to each other
    def perp2D(u):
        return array([u[1], -u[0]])

    # check if the S is just one single point
    if ((S[0][0] == S[1][0]) and (S[0][1] == S[1][1])):
        cn = cn_PnPoly(S[0], V)
        if (cn == 1):
            return S
        else:
            return None

    tE = 0.
    tL = 1.

    dS = S[1] - S[0]

    n = len(V)

    for i in range(n):
        e = V[(i+1)%n] - V[i]
        ne = perp2D(e)
        N = -dot(ne, S[0]-V[i])
        D = dot(ne, dS)

        if (fcmp(D, 0)):
            if (N < 0):
                return None
            else:
                continue

        t = N / D
        if (D < 0):
            if (t > tE):
                tE = t
                if (tE > tL):
                    return None
        else:
            if (t < tL):
                tL = t
                if (tL < tE):
                    return None

    int_s = [S[0] + tE * dS, S[0] + tL*dS]
        
    return int_s

#
# Function intersectPolyPoly2D
#
# Find the intersection between two 2D polygons
#
def intersectPolyPoly2D(p1, p2):
    segs = []

    n1 = len(p1)
    n2 = len(p2)

    # if the two polygons are equal, i.e. have the same number of vertices and
    # all vertices are shared, then simply return that polygon
    if (n1 == n2):
        a = 0
        for i in range(n1):
            if not ((fcmp(p1[i][0], p2[i][0])) and
                    (fcmp(p1[i][1], p2[i][1]))):
                a = 1
                break
        if (a == 0):
            return p1

    # now we compute the intersections between the segments of polygon 1
    # and polygon 2 and vice versa
    # this way no vertex of the intersection is missed, but each point appears
    # twice, in the end of segment s(i) and also in the beginning of s(i+1)
    for i in range(n1):
        s = [p1[i], p1[(i+1)%n1]]
        iseg = intersectSegPoly2D(s, p2)
        
        # if none then the line segment doesn't intersect with the polygon
        if (iseg is not None):
            # make sure that the line segment is not a single point
            if not ((fcmp(iseg[0][0], iseg[1][0])) and
                    (fcmp(iseg[0][1], iseg[1][1]))):
                segs.append(iseg)

    for i in range(n2):
        s = [p2[i], p2[(i+1)%n2]]
        iseg = intersectSegPoly2D(s, p1)

        # if none then the line segment doesn't intersect with the polygon
        if (iseg is not None):
            # make sure that the line segment is not a single point
            if not ((fcmp(iseg[0][0], iseg[1][0])) and
                    (fcmp(iseg[0][1], iseg[1][1]))):
                segs.append(iseg)

    # this takes care of the double counting of vertices
    if (len(segs) > 0):
        int_poly = polygonFromSeg2D(segs)
        return int_poly

    return None


#
# Function intersectSegPolyhedron3D
#
# Find the intersections between a line segment s (two points) and a convex 3D
# polyhedron {V,F} (n vertices, n faces, where a face is a normal vector)
# (taken from: http://softsurfer.com/Archive/algorithm_0111/algorithm_0111.htm)
# Returns either the segment that intersects or None if error
#
def intersectSegPolyhedron3D(S, V, F):
    # check if the s is just one single point
    if (fcmp(S[0][0], S[1][0]) and fcmp(S[0][1], S[1][1]) and fcmp(S[0][2], S[1][2])):
        cn = cn_PnPolyhedron(S[0], V, F)
        if (cn == 1):
            return S
        else:        
            return None

    tE = 0.
    tL = 1.
    
    dS = S[1] - S[0]

    n = len(V)

    for i in range(n):
        N = -dot(F[i], S[0]-V[i])
        D = dot(F[i], dS)

        if (fcmp(D, 0)):
            if (N < 0):
                return None
            else:
                continue
            
        t = N / D
        if (D < 0):
            if (t > tE):
                tE = t
                if (tE > tL):
                    return None
        else:
            if (t < tL):
                tL = t
                if (tL < tE):
                    return None

    int_s = [S[0] + tE*dS, S[0] + tL*dS]

    return int_s

#
# Function intersectPolyhedra3D
#
# Find the intersection between two polyhedra
# each polyhedron is given as a set of line segments S, and a set of faces, each face consists
# of a set of vertices and a set of normal vectors (equally sized)
# Returns a set of line segments, which points define the interecting volume
#
def intersectPolyhedra3D(S1, V1, F1, S2, V2, F2):
    segs = []

    n1 = len(S1)
    n2 = len(S2)

    # find all intersections between segments of polyhedron 1 (S1) and polyhedron 2 (V2, F2)
    # and vice versa
    # this results in a set of line segments that represent the interected volume, though
    # vertices may have multiplicity 
    for i in range(n1):
        iseg = intersectSegPolyhedron3D(S1[i], V2, F2)
        if (iseg is not None):
            if not (fcmp(iseg[0][0], iseg[1][0]) and
                    fcmp(iseg[0][1], iseg[1][1]) and
                    fcmp(iseg[0][2], iseg[1][2])):
                segs.append(iseg)

    for i in range(n2):
        iseg = intersectSegPolyhedron3D(S2[i], V1, F1)
        if (iseg is not None):
            if not (fcmp(iseg[0][0], iseg[1][0]) and
                    fcmp(iseg[0][1], iseg[1][1]) and
                    fcmp(iseg[0][2], iseg[1][2])):
                segs.append(iseg)

    if (len(segs) > 0):
        int = polyhedronCompact3D(segs)
        return int
    
    return None

# Test 1: the 3D polyhedron intersection function
def test1_intPolyhedra():
    # vertices of polyhedron p
    p1 = array([1., 0., 0.])
    p2 = array([1., 1., 0.])
    p3 = array([0., 1., 0.])
    p4 = array([0., 0., 0.])
    p5 = array([1., 0., 1.])
    p6 = array([1., 1., 1.])
    p7 = array([0., 1., 1.])
    p8 = array([0., 0., 1.])

    # vertices of polyhedron q
    q1 = array([1., 0., 0.])
    q2 = array([1., 1., 0.])
    q3 = array([0., 1., 0.])
    q4 = array([0., 0., 0.])
    q5 = array([1., 0., 1.])
    q6 = array([1., 1., 1.])
    q7 = array([0., 1., 1.])
    q8 = array([0., 0., 1.])

    # normal vectors
    n1 = array([1., 0., 0.])
    n2 = array([-1., 0., 0.])
    n3 = array([0., 1., 0.])
    n4 = array([0., -1., 0.])
    n5 = array([0., 0., 1.])
    n6 = array([0., 0., -1.])

    # polyhedron 1: sides (S1), faces (V1,F1)
    S1 = [[p1,p2], [p2,p3], [p3,p4], [p4,p1], [p1,p5], [p2,p6], [p3,p7], [p4,p8],
          [p5,p6], [p6,p7], [p7,p8], [p8,p5]]
    V1 = [p6, p6, p6, p4, p4, p4]
    F1 = [n1, n3, n5, n2, n4, n6]

    # polyhedron 2: sides (S2), faces (V2,F2)
    S2 = [[q1,q2], [q2,q3], [q3,q4], [q4,q1],
          [q1,q5], [q2,q6], [q3,q7], [q4,q8],
          [q5,q6], [q6,q7], [q7,q8], [q8,q5]]
    V2 = [q6, q6, q6, q4, q4, q4]
    F2 = [n1, n3, n5, n2, n4, n6]

    v = intersectPolyhedra3D(S1, V1, F1, S2, V2, F2)
    print len(v)
    for i in range(len(v)):
        print v[i]

    print polyhedronVolume3D(v)

# Test 2: the 2D polygon intersection function
def test2_intPoly():
    p1 = array([1., 0.])
    p2 = array([1., 1.])
    p3 = array([0., 1.])
    p4 = array([0., 0.])    

    q1 = array([1., 0.5])
    q2 = array([1., 1.])
    q3 = array([0.5, 1.])
    q4 = array([0.5, 0.5])

    P1 = [p1, p2, p3, p4]
    P2 = [q1, q2, q3, q4]

    segs = intersectPolyPoly2D(P1, P2)
    print len(segs)
    print segs

def main(argv = None):
    if (argv is None):
        argv = sys.argv
    argc = len(argv)
    
    if (argc > 1):
        test = int(argv[1])
    else:
        test = 0
    
    if (test == 1):
        # Test 1: intersection between two polyhedrons
        print "Test 1: intersection between two polyhedrons"
        test1_intPolyhedra()
    elif (test == 2):
        # Test 2: intersection between two polygons
        print "Test 2: intersection between two polygons"
        test2_intPoly()
    elif (test == 3):
        # Test 3: test the cartesian to GC converter
        print "Test 3: test the cartesian to GC converter"
        for phi in range(-18,19):
            x,y,z = convertGCToCartesian(phi*10., 0.)
            l,b = convertCartesianToGC(x, y, z)
            print l,b
    elif (test == 4):
        # Test 4: test the rotation by Euler angles
        print "Test 4: test the rotation by Euler angles"
        
        l = 120.0
        b = -30.0
        u = array(convertGCToCartesian(l, b))
        print u

        x = array([1.0, 0.0, 0.0])
        y = array([0.0, 1.0, 0.0])
        z = array([0.0, 0.0, 1.0])

        phi = -90.0
        theta = -b
        psi = 90.0 - l
        print "phi =", phi, "; theta =", theta, "; psi =", psi
        xp = array(rotate3D(x, rad(phi), rad(theta), rad(psi)))
        print "x' =", xp
        print "x =", array(rotate3Dinverse(xp, rad(phi), rad(theta), rad(psi)))
        yp = array(rotate3D(y, rad(phi), rad(theta), rad(psi)))
        print "y' =", yp
        zp = array(rotate3D(z, rad(phi), rad(theta), rad(psi)))
        print "z' =", zp

if __name__ == "__main__":
    sys.exit(main())
