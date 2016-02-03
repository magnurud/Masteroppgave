#! /bin/bash
from numpy import *
from writepoints import *

# A script that generates the hpts.in file to be used in nek. 
# generates the wanted points equidistantly distributed along a line

rho = 100; # Number of lines used to create the planes.
np = 200   # Number of points in each line
nl = 12    # number of lines

N = np*nl+2*rho**2 # Total number of history points

# x is the array containing the starting and ending values of the line
# x = x0,x1,y0,y1,z0,z1

x = array([[-2.62,-2.62, -1.00,1.00,0.025,0.025],#horizontal line 
            [0.000,0.000,-1.00,1.00,0.025,0.025],#horizontal line 
            [2.420,2.420,-1.00,1.00,0.025,0.025],#horizontal line 
            [4.260,4.260,-1.00,1.00,0.025,0.025],#horizontal line 
            [8.840,8.840,-1.00,1.00,0.025,0.025],#horizontal line 
            [18.02,18.02,-1.00,1.00,0.025,0.025],#horizontal line 
            [-2.62,-2.62,0.000,0.00,0.000,0.500],#vertical line 
            [0.000,0.000,0.000,0.00,0.000,0.500],#vertical line 
            [2.420,2.420,0.000,0.00,0.000,0.500],#vertical line 
            [4.260,4.260,0.000,0.00,0.000,0.500],#vertical line 
            [8.840,8.840,0.000,0.00,0.000,0.500],#vertical line 
            [18.02,18.02,0.000,0.00,0.000,0.500],#vertical line 
            ])
x[:,[0,1]] = x[:,[0,1]]*0.109+1.4315 # Scaling the xcoords correctly
        
writenumberofpoints(N) # initialize hpts.in with number of points

s = 'a' # From now on we will only append to hpts.in

# Writing the coordinates for the lines to hpts.in 
for i in range(nl):
	makehpts2(genline,np,x[i],s)
	#makehpts(genline,np,x[i],s)

# Writing the coordinates for the full planes to hpts.in 
for zval in linspace(0.000,1.5,rho+2):
    if(zval != 0.000 and zval!= 1.5):
        x = array([0.001,0.001, -1.00,1.00,zval,zval]) #horizontal lines
        makehpts2(genline,rho,x,s)

for zval in linspace(0.000,1.5,rho+2):
    if(zval != 0.000 and zval!= 1.5):
        x = array([0.300,0.300, -1.00,1.00,zval,zval]) #horizontal lines
        makehpts2(genline,rho,x,s)
