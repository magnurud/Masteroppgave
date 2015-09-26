#! /bin/bash
from numpy import *

def points2circ(x1,x2,x3):
    # A function that takes in three 3D points and returns 
    # The radius (r) and center x of the sphere with its center 
    # in the same plane as the three points.
    # INPUT: xi = array([real,real,real])
    # OUTPUT: r = real, x = array([real,real,real])
    n = cross(x2-x1,x3-x1) # 
    phi1 = x2-x1
    psi1 = -sum(x1**2-x2**2)/2
    phi2 = x3-x2
    psi2 = -sum(x2**2-x3**2)/2
    psi3 = dot(n,x1)
    A = array([phi1,phi2,n])
    b = array([psi1,psi2,psi3])
    x = linalg.solve(A,b)
    r = sqrt(dot(x-x1,x-x1))
    return r,x
