#! /bin/bash
from sympy import *
import numpy as np

# Different Constants
a = Symbol("a")    # Constants for different uses
b = Symbol("b")    # Constants for different uses
c = Symbol("c")    # Constants for different uses
d = Symbol("d")    # Constants for different uses
e = Symbol("e")    # Constants for different uses
# Corner values 
c1= Symbol("c1")    # Constants for different uses
c2= Symbol("c2")    # Constants for different uses
c3= Symbol("c3")    # Constants for different uses
c4= Symbol("c4")    # Constants for different uses
# Reference element values 
xi = Symbol("xi")    # Xi 
eta = Symbol("eta")    # Eta
t = Symbol("t") # Theta angle 
# Basis functions
p = Symbol("p") # basis functions for xi 
q = Symbol("q") # basis functions for xi 
r = Symbol("r") # basis functions for et
# Gordon hall functions
Fe = Symbol("Fe") # First contribution
Fn = Symbol("Fn") # First contribution
Fen = Symbol("Fen") # First contribution
Fgh = Symbol("Fgh") # First contribution

# Defining edge functions 
gamma1 = Matrix([a +c*(1+eta)/2,b*c*(1+eta)/2])
gamma2 = Matrix([-xi*(a+c),b*c])# + k*(1-xi**2)
gamma3 = Matrix([-a-c*(1+eta)/2,b*c*(1+eta)/2])
gamma4 = Matrix([cos(t), sin(t)-sin(np.pi/2-e)])

#### HAVE TO CHECK THIS !!! 
# Evaluating in the corners
expr1 = gamma1; expr2 = gamma2; expr3 = gamma3; expr4 = gamma4;
c1 = expr1.subs(eta,-1)
c2 = expr2.subs(eta,-1)
c3 = expr3.subs(eta,-1)
c4 = expr4.subs(eta,-1)

# Defining basis functions 
p = (1-xi)/2
q = (1+xi)/2
r = (1-eta)/2
s = (1+eta)/2

#Defining the Mapping 
Fe = p*gamma3+q*gamma1
Fn = r*gamma4+s*gamma2
Fen = p*r*c1 + p*s*c2 + q*r*c3 + q*s*c4

Fgh = Fe+Fn-Fen

#pprint(Fe)
#pprint(Fn)
#pprint(Fgh[0])
# Calculating Jacobi-Matrix
J = Matrix([[diff(Fgh[0],xi),diff(Fgh[0],eta)],
    [diff(Fgh[1],xi),diff(Fgh[1],eta)]])

# Caldulating determinant and changes in determinant
Jdet = J[0,0]*J[1,1]-J[0,1]*J[1,0]
Jdiv1 = diff(Jdet,xi) 
Jdiv2 = diff(Jdet,eta) 
