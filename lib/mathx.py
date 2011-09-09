#
# Python/numarray repixeling package
#   module: mathx.py
#

from numarray import *

# machine epsilon used for floating point comparisons,
# i.e. (x - y < EPSILON) instead of (x == y)
EPSILON = 1e-15

# --------------------------------------------------------------------------------
#  Function definitions
# --------------------------------------------------------------------------------

#
# Function fcmp
#
# Compares floats
#
def fcmp(x, y):
    d = fabs(x - y)
    if (d < EPSILON):
        return 1
    else:
        return 0

#
# Function nacmp
#
# Compares numarrays
#
def nacmp(a, b):
    if (a.getshape() != b.getshape()):
        return 0

#
# Function rad
#
# Convert degrees to radians
#
def rad(x):
    return x * pi/180.

#
# Function deg
#
# Convert radians to degrees
#
def deg(x):
    return x * 180./pi

#
# Function crossprod
#
# Calculate cross product between two vectors
#
def crossprod(u, v):
    n = array([u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]])
    return n
