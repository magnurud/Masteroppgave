#! /bin/bash
from numpy import *

def genline(t,x):
    # INPUT: 
    # t - parameter
    # x contains 6 values
        # x0 - starting point for x-value
        # x1 - ending point for x-value
        # y0 - starting point for y-value
        # y1 - ending point for y-value
        # z0 - starting point for z-value
        # z1 - ending point for z-value
    # OUTPUT:
    # a function mapping from 0,1 to a set of 3D points
    return array([[t*(x[1]-x[0])+x[0]],[t*(x[3]-x[2])+x[2]],[t*(x[5]-x[4])+x[4]]])

def makehpts(f,n,x,s):
    # INPUT: 
    # f - function which defines a line
    # n - number of points
    # s - a string if s==a => append , if s==w => rewrite
    # XYZ - an array of length 6 that contains the endpoint coordinates
    # OUTPUT:
    # a file hpts.in containing n the points on the line f(x), x in (0,1)
    with file('hpts.in',s) as outfile:
        for t in linspace(0,1,n):
         #   outfile.write(f(t,x))
            savetxt(outfile,transpose(f(t,x)),fmt='%-7.2f')
        #outfile.close()
    return 1

def makehpts2(f,n,x,s):
    # FOR INTERNAL POINTS ONLY !! 
    # INPUT: 
    # f - function which defines a line
    # n - number of points
    # s - a string if s==a => append , if s==w => rewrite
    # XYZ - an array of length 6 that contains the endpoint coordinates
    # OUTPUT:
    # a file hpts.in containing n the points on the line f(x), x in (0,1)
    with file('hpts.in',s) as outfile:
        for t in linspace(0,1,n+2):
         #   outfile.write(f(t,x))
         if(t != 0 and t!= 1):
            savetxt(outfile,transpose(f(t,x)),fmt='%-7.3f')
        #outfile.close()
    return 1
