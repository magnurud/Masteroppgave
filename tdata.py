# /bin/bash
from numpy import *

#Nek performance
x = array([
8.44357,
5.31558E+00,
2.77470E+00,
1.95860E+00,
1.67716E+00
])
n = array([1,2,4,6,8])

#Fluent performance
x = array([
18.9722, #20
10.804, #40
3.4583, #100
2.9266  #120
])
n = array([1,2,5,6])

#Fluent performance per iteration
x = array([
18.9722,#/(50-5.6740), #20
10.804, #/(50-5.6298), #40
3.4583,#/(50-14.1771), #100
2.9266#/(50-14.0845)  #120
])
n = array([1,2,5,6])

x = x/x[0]#*n

print 1/x


#[ 1.          0.95786572  0.90465639  0.79234837]
