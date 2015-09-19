import numpy as np
from element import plotjacobi

case = 3
angle =np.pi/8  # half the size of the total angle 
c1 = np.cos(np.pi/2-2*angle) # lower Corner x-values
c2 = np.tan(np.pi/2-2*angle) # Slope for the linear boundaries
c3 = 0.3 # lateral expansion
k = 0.7*(1-np.sin(np.pi/2-angle)) # 2nd order coefficient

if case <= 3:
    def gamma1(eta):
        return c1+c3*(1+eta)/2,c2*c3*(1+eta)/2
        #return 1,eta # Linear
    if case == 1:
        def gamma2(xi):
            return -xi*(c1+c3),c2*c3# + k*(1-xi**2)
            #return xi,1 # Linear
    elif case == 2:
			def gamma2(xi):
				theta = -angle*xi+np.pi/2     
				return c3**2*(1+c2**2)*np.cos(theta),c3**2*(1+c2**2)*np.sin(theta)-np.sin(np.pi/2-angle)
            #return xi,1 # Linear
    elif case == 3:
        def gamma2(xi):
            return -xi*(c1+c3),c2*c3 + k*(1-xi**2)
            #return xi,1 # Linear
    def gamma3(eta):
        return -c1-c3*(1+eta)/2,c2*c3*(1+eta)/2
        #return -1,eta # Linear

    def gamma4(xi):
        theta = -angle*xi+np.pi/2     
        return c1*xi, np.sin(theta)-np.sin(np.pi/2-angle)
        #return np.cos(theta), np.sin(theta)-np.sin(np.pi/2-angle)
        #return xi,-1 # Linear
if case == 4 or case == 5:
    def gamma1(eta):
        return c1,c2*c3*(1+eta)/2
        #return 1,eta # Linear
    if case == 4:
        def gamma2(xi):
            return -xi*c1,c2*c3# + k*(1-xi**2)
            #return xi,1 # Linear
    elif case == 5:
        def gamma2(xi):
            return -xi*c1,c2*c3+k*(1-xi**2)
            #return xi,1 # Linear
    def gamma3(eta):
        return -c1,c2*c3*(1+eta)/2
        #return -1,eta # Linear

    def gamma4(xi):
        theta = -angle*xi+np.pi/2     
        return np.cos(theta), np.sin(theta)-np.sin(np.pi/2-angle)
        #return xi,-1 # Linear
plotjacobi(gamma1,gamma2,gamma3,gamma4)
