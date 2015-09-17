#! /bin/bash
from numpy import *
from writepoints import *
# A script that generates the hpts.in file to be used in nek. 
# generates the wanted points equidistantly distributed along a line
np = 10 # Number of points in each line
nl = 12 # number of lines
# x is the array containing the starting and ending values of the line
# x = x0,x1,y0,y1,z0,z1

x = array([[-2.62,-2.62,-1.75,1.75,0.025,0.025], #horizontal lines
					 [0.000,0.000,-1.75,1.75,0.025,0.025],#horizontal lines
					 [2.420,0.420,-1.75,1.75,0.025,0.025],#horizontal lines
					 [4.260,4.260,-1.75,1.75,0.025,0.025],#horizontal lines
					 [8.840,8.840,-1.75,1.75,0.025,0.025],#horizontal lines
					 [18.02,18.02,-1.75,1.75,0.025,0.025],#horizontal lines
					 [-2.62,-2.62,0.000,0.00,0.000,1.500], #vertical lines
					 [0.000,0.000,0.000,0.00,0.000,1.500],#vertical lines
					 [2.420,2.420,0.000,0.00,0.000,1.500],#vertical lines
					 [4.260,4.260,0.000,0.00,0.000,1.500],#vertical lines
					 [8.840,8.840,0.000,0.00,0.000,1.500],#vertical lines
					 [18.02,18.02,0.000,0.00,0.000,1.500],#vertical lines
					])
s = 'w'
for i in range(nl):
  # makehpts2 FOR INTERNAL POINTS ONLY !! 
	makehpts2(genline,np,x[i],s)
	s = 'a'

