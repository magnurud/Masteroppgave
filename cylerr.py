#! /bin/bash
from numpy import *

cl = 0.009401
cd = 6.18533

x = array([[ 6.1845124E+00, 7.0473316E-03 ], #horizontal lines
        [6.1853362E+00 ,8.2500749E-03 ],
        [6.1852302E+00 ,8.2500196E-03],
        [6.1331543E+00 ,8.2137267E-03],
        [6.1795558E+00, 8.8143970E-03],
        [6.2101968E+00, 6.4789122E-03], # Default p = 6
        [6.2104179E+00, 6.4658555E-03] # ifchar = T


        ])

#print x[1,:]
print x

x[:,0] = x[:,0]-cd 
x[:,0] = x[:,0]/cd*100 

x[:,1] = x[:,1]-cl 
x[:,1] = x[:,1]/cl*100 

print x
