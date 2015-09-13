import numpy as np
import matplotlib.pyplot as plt
import gh
import quadrature_nodes as qn
import structured_grids as sg 
import laplace_functions as lp 


#N = 10
#xis = np.linspace(-1, 1, num=N)

#Number of GLL points:
N = 20
xis = qn.GLL_points(N)
weights = qn.GLL_weights(N, xis)
etas = xis
k = 0.5
angle =3*np.pi/8

def gamma1(eta):
    #theta = np.pi/4.*(eta+1)     
    #return 2.*np.cos(theta), 2.*np.sin(theta)
    #return (1-eta),1+eta
    return 1/np.sqrt(2)+(eta+1)/2,np.tan(angle)*(1+eta)/2


def gamma2(xi):
    #return -xi*(1/np.sqrt(2)+1),np.tan(angle)*(1+xi)/2 + 0.5*(1-xi**2)
    return -xi*(1/np.sqrt(2)+1),np.tan(angle)*(1+xi)/2# + 0.5*(1-xi**2)

def gamma3(eta):
    #theta = np.pi/8.*(eta)
    #return np.cos(theta), np.sin(theta)
    return -1/np.sqrt(2)+(-eta+1)/2,np.tan(angle)*(1+eta)/2

def gamma4(xi):
    theta = (np.pi/2-angle)*(xi)+np.pi/2     
    return np.cos(theta), np.sin(theta)
    #return 1.0/np.sqrt(2)*xi,0

#Generate mesh:
print "Getting Mesh..."
X,Y = gh.gordon_hall_grid(gamma1, gamma2, gamma3, gamma4, xis, xis)

# Get lagrangian derivative matrix:
print "Getting derivative matrix..."
D = sg.diff_matrix(xis, N)

# Get Jacobian and total G-matrix:
print "Getting Jacobian and G-matrices."
Jac, G_tot = lp.assemble_total_G_matrix(np.dot(X,D),
                                        np.dot(X.T,D),
                                        np.dot(Y,D),
                                        np.dot(Y.T,D),
                                        N,
                                        N)
#
print "Assembling stiffness matrix."
A = lp.assemble_local_stiffness_matrix(D, G_tot, N, weights)

#print Jac
print np.linalg.cond(A)
#print G_tot

#Plot all rows:
#for i in range(N):
    #plt.plot(X[i,:],Y[i,:],'b')
    #plt.plot(X[:,i],Y[:,i],'b')
#plt.show()
