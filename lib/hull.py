#
# Python/numarray repixeling package
#   module: hull.py
#

from numarray import *
from mathx import *

# --------------------------------------------------------------------------------
#  Class definitions
# --------------------------------------------------------------------------------

# We do not define classes for points, vertices and vectors. Plain numarrays
# are used for these objects

#
# Class Edge
#
# Representation of an edge on a polytope/polyhedron (line between two vertices)
#
class Edge:

    #
    # Constructor
    #
    def __init__(self, u, v, facet):
        # u, v are the two vertices

        self.u = u
        self.v = v
        self.facet = facet

        self.next = None
        self.prev = None
        self.twin = None

    #
    # Method matches
    #
    def matches(self, u, v):
        return ((fcmp(self.u[0], u[0]) and fcmp(self.u[1], u[1]) and fcmp(self.u[2], u[2]) and
                 fcmp(self.v[0], v[0]) and fcmp(self.v[1], v[1]) and fcmp(self.v[2], v[2])) or
                (fcmp(self.u[0], v[0]) and fcmp(self.u[1], v[1]) and fcmp(self.u[2], v[2]) and
                 fcmp(self.v[0], u[0]) and fcmp(self.v[1], u[1]) and fcmp(self.v[2], u[2])))

    #
    # Method onHorizon
    #
    def onHorizon(self):
        if (self.twin is None):
            return 0
        else:
            return (self.facet.marked == 0) and (self.twin.facet.marked == 1)

    #
    # Method findHorizon
    #
    def findHorizon(self, horizon):
        #print "--- findHorizon ---"
        if (self.onHorizon()):
            if (len(horizon) > 0) and (self == horizon[0]):
                #print "--- end findHorizon ---"
                return
            else:
                horizon.append(self)
                self.next.findHorizon(horizon)
        else:
            if (self.twin is not None):
                self.twin.next.findHorizon(horizon)

        #print "--- end findHorizon ---"
                
#
# Class Facet
#
# Representation of a facet on a polytope/polyhedron (always in 3D)
#
class Facet:

    #
    # Constructor
    #
    def __init__(self, v, ref = None):
        # v is a list of coplanar vertices
        # ref is a point used for orientation

        self.v = v
        self.marked = 0    # flag mark if to be deleted
        n = len(v)
        self.normal = crossprod(self.v[1]-self.v[0], self.v[n-1]-self.v[0])

        if (ref is not None):
            self.orient(ref)

        self.createEdges()

    #
    # Method behind
    #
    def behind(self, ref):
        return (dot(self.normal, ref) < dot(self.normal, self.v[0]))

    #
    # Method connectFacet
    #
    def connectFacet(self, adjacent, u, v):
        inner = self.getMatchingEdge(u, v)
        outer = adjacent.getMatchingEdge(u, v)

        inner.twin = outer
        outer.twin = inner

    #
    # Method connectEdge
    #
    def connectEdge(self, e):
        inner = self.getMatchingEdge(e.u, e.v)
        inner.twin = e
        e.twin = inner

    #
    # Method getMatchingEdge
    #
    def getMatchingEdge(self, u, v):
        n = len(self.e)
        for i in range(n):
            if (self.e[i].matches(u, v)):
                return self.e[i]
            
        return None

    #
    # Method getHorizonEdge
    #
    def getHorizonEdge(self):
        #print "--- getHorizonEdge ---"
        n = len(self.e)
        for i in range(n):
            opposite = self.e[i].twin
            #print self.e[i].u, self.e[i].v
            #print opposite.u, opposite.v
            if (opposite is not None):
                if (opposite.onHorizon()):
                    #print "-+- end getHorizonEdge -+-"
                    return self.e[i]
                
        #print "--- end getHorizonEdge ---"
        return None

    #
    # Method orient
    #
    def orient(self, ref):
        if (not self.behind(ref)):
            tmp = [self.v[0]]
            n = len(self.v)
            for i in range(n-1):
                tmp.append(self.v[n-i-1])
            del self.v
            self.v = tmp
            
            self.normal = -self.normal
            self.createEdges()

    #
    # Method createEdges
    #
    def createEdges(self):
        self.e = []
        n = len(self.v)
        for i in range(n):
            edge = Edge(self.v[i], self.v[(i+1)%n], self)
            self.e.append(edge)

        n = len(self.e)
        for i in range(n):
            self.e[i].next = self.e[(i+1)%n]
            self.e[i].prev = self.e[(i-1)%n]

#
# Class ConvexHull
#
# Note: When the first four points are co-planar, this algorithm fails!
#       How fix this?
#
class ConvexHull3D:

    #
    # Constructor
    #
    def __init__(self):
        self.facets = []
        
        self.created = []
        self.horizon = []
        self.visible = []

    #
    # Method create
    #
    # Create the convex hull
    #
    def create(self, v):
        self.prepare(v)
        n = len(v)
        ##print "remaining points: ", n-4
        for i in range (n-4):
            self.stepA(v)
            self.stepB(v)
            self.stepC(v)

    #
    # Method prepare
    #
    def prepare(self, v):
        # start out by creating a tetrahedron
        
        # we must make sure that the four points chosen does in fact
        # form a tetrahedron!!!

        v1 = v[1] - v[0]
        v2 = v[2] - v[0]
        n = crossprod(v1, v2)
        ip = 0
        i = 4
        while (ip == 0):
            v3 = v[3] - v[0]
            # calculate the dot product to see if on the same plane
            # if not, break, otherwise interchange two vertices
            if (innerproduct(v3, n) != 0):
                break
            else:
                tmp = v[3]
                v[3] = v[i]
                v[i] = tmp
                i += 1
        
        f1 = Facet([v[0], v[1], v[2]], ref = v[3])
        self.facets.append(f1)
        
        f2 = Facet([v[0], v[2], v[3]], ref = v[1])
        self.facets.append(f2)
        
        f3 = Facet([v[0], v[1], v[3]], ref = v[2])
        self.facets.append(f3)
        
        f4 = Facet([v[1], v[2], v[3]], ref = v[0])
        self.facets.append(f4)

        f1.connectFacet(f2, v[0], v[2])
        f1.connectFacet(f3, v[0], v[1])
        f1.connectFacet(f4, v[1], v[2])
        f2.connectFacet(f3, v[0], v[3])
        f2.connectFacet(f4, v[2], v[3])
        f3.connectFacet(f4, v[1], v[3])

        self.current = 4

    #
    # Method stepA
    #
    def stepA(self, v):
        #print "--- stepA ---"
        if (self.current >= len(v)):
            return

        self.created = []
        self.horizon = []
        self.visible = []

        #print "current: ", v[self.current]
    
        n = len(self.facets)
        for i in range(n):
            f = self.facets[i]
            #print "facet       : ", f.v
            #print "facet normal: ", f.normal
            #print "      behind: ", f.behind(v[self.current])
            if (f.behind(v[self.current]) == 0):
                self.visible.append(f)

        # if point is inside the current hull then try next point
        if (len(self.visible) == 0):
            self.current += 1
            self.stepA(v)
            return

        #print "no visible: ", len(self.visible)

        # mark all visible facets for deletion
        n = len(self.visible)
        for i in range(n):
            self.visible[i].marked = 1

        # find horizon edges
        for i in range(n):
            e = self.visible[i].getHorizonEdge()
            if (e is not None):
                #print "horizon edge: ", e.u, e.v
                #print "        twin: ", e.twin.u, e.twin.v
                #print "   onHorizon: ", e.onHorizon()
                #print "    f marked: ", e.facet.marked
                #print "  t f marked: ", e.twin.facet.marked
                #print "(b) horizon: ", self.horizon
                e.findHorizon(self.horizon)
                #print "(a) horizon: ", self.horizon
                break
        
    #
    # Method stepB
    #
    def stepB(self, v):
        #print "--- stepB ---"
        if (self.current >= len(v)):
            return

        old = None
        last = None
        first = None

        n = len(self.horizon)
        #print " no horizon: ", n
        for i in range(n):
            e = self.horizon[i]
            old = e.twin.facet

            # create new facet
            f = Facet([v[self.current], e.v, e.u])
            self.facets.append(f)
            self.created.append(f)

            f.connectEdge(e)
            if (last is not None):
                f.connectFacet(last, v[self.current], e.u)
            last = f
            if (first is None):
                first = f

        if ((last is not None) and (first is not None)):
            last.connectFacet(first, v[self.current], first.e[1].v)

        #print "created: ", self.created
        
    #
    # Method stepC
    #
    def stepC(self,v ):
        #print "--- stepC ---"
        n = len(self.visible)
        for i in range(n):
            f = self.visible[i]
            self.facets.remove(f)

        self.created = []
        self.current += 1

    #
    # Method getFacets
    #
    def getFacets(self):
        list = []
        n = len(self.facets)
        for i in range(n):
            list.append(self.facets[i].v)

        return list

# --------------------------------------------------------------------------------
#  Function definitions
# --------------------------------------------------------------------------------

#
# Function convexhull_3D
#
# Constructs a 3D convex hull from a set of n points in R3
# and return a list of triangular facets
#
def convexHull3D(V):
    hull = ConvexHull3D()
    hull.create(V)
    return hull.getFacets()
