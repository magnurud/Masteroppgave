#! /bin/bash
from sympy import *
import numpy as np

#    _<__(2)___
#   |          |
#   |          ^
#   |          |
#  (3)        (1)
#   |          |
#   |___(4)_>__| 
#
#
# Different Constants
a = Symbol("a")    # Constants for different uses
b = Symbol("b")    # Constants for different uses
c = Symbol("c")    # Constants for different uses
d = Symbol("d")    # Constants for different uses
e = Symbol("e")    # Constants for different uses
# Corner values 
c1,c2,c3,c4 = symbols('c1,c2,c3,c4')    # Constants for corner values
# Reference element values 
xi,eta,theta = symbols('xi,eta,theta')    # Xi 
# Basis functions
p,q,r,s = symbols('p,q,r,s') # basis functions for xi and eta
# Gordon hall functions
Fe,Fn,Fen,Fgh = symbols("Fe,Fn,Fen,Fgh") # Gordon Hall parts

angle = pi/4  # half the size of the total angle 
a = cos(pi/2-angle) # lower Corner x-values
b = tan(pi/2-angle) # Slope for the linear boundaries
c = 0.3 # lateral expansion
d = 1+c # lateral endpoints
#d = sin(np.pi/2-angle)
#k = 0.7*(1-np.sin(np.pi/2-angle)) # 2nd order coefficient
k = 0.3
theta = pi/4*(2-xi)    

# Defining edge functions 
### trapezoidal element with curved bottom
#gamma1 = Matrix([a +c*(1+eta)/2,b*c*(1+eta)/2])
##gamma2 = Matrix([-xi*(a+c),b*c])# + k*(1-xi**2)
#gamma2 = Matrix([d*cos(theta),d*sin(theta)-sin(pi/4)])
#gamma3 = Matrix([-a-c*(1+eta)/2,b*c*(1+eta)/2])
#gamma4 = Matrix([cos(theta),sin(theta)-sin(pi/4)])
#gamma4 = Matrix([a*xi,sqrt(1-(a*xi)**2)])
### Circular element 
#a = 1/np.sqrt(2) # cornervals
#gamma1 = Matrix([sqrt(1-(a*eta)**2),a*eta])
#gamma2 = Matrix([-a*xi,sqrt(1-(a*xi)**2)])
#gamma3 = Matrix([-sqrt(1-(a*eta)**2),-a*eta])
#gamma4 = Matrix([a*xi,-sqrt(1-(a*xi)**2)])
### quadratic element with curved bottom (and top
gamma1 = Matrix([1,eta])
gamma2 = Matrix([-xi,1 + 0*k/2*(1-xi**2)])
#gamma2 = Matrix([-xi,1])
gamma3 = Matrix([-1,-eta])
gamma4 = Matrix([xi,-1+0*k*(1-xi**2)])
#gamma4 = Matrix([a*xi,sqrt(1-(a*xi)**2)])


#### HAVE TO CHECK THIS !!! 
# Evaluating in the corners
expr1 = gamma1; expr2 = gamma2; expr3 = gamma3; expr4 = gamma4;
c1 = expr4.subs(xi,-1); c2 = expr3.subs(eta,-1) # left side corners
c3 = expr1.subs(eta,-1); c4 = expr2.subs(xi,-1) # right side corners
# Defining basis functions 
p = (1-xi)/2; q = (1+xi)/2
r = (1-eta)/2; s = (1+eta)/2

#Defining the Mapping 
Fe = p*gamma3+q*gamma1
Fn = r*gamma4+s*gamma2
Fen = p*r*c1 + p*s*c2 + q*r*c3 + q*s*c4

Fgh = Fe+Fn-Fen

#pprint(Fe)
#pprint(Fn)
#pprint(simplify(Fgh[0]))
# Calculating Jacobi-Matrix
J = Matrix([[diff(Fgh[0],xi),diff(Fgh[0],eta)],
    [diff(Fgh[1],xi),diff(Fgh[1],eta)]])

pprint(J)
# Caldulating determinant and changes in determinant
Jdet = simplify(J[0,0]*J[1,1]-J[0,1]*J[1,0])
Jdiv1 = diff(Jdet,xi) 
#Jdiv2 = diff(Jdet,eta) 
#pprint(collect(Jdet,eta))
pprint(Jdet)

