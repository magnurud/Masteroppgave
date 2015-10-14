from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import gh
import quadrature_nodes as qn
import structured_grids as sg 
import laplace_functions as lp 
#Number of GLL points:
N = 11
xis = qn.GLL_points(N)
one = np.ones(N)
weights = qn.GLL_weights(N, xis)
etas = xis
pi = 4*np.arctan(1)
#Plotting the gll-mesh
#for i in range(N):
	#plt.plot(xis[i]*one,xis[:],'b')
	#plt.plot(xis[:],xis[i]*one,'b')

# Circle
M = 100
circle = np.zeros((M,2))
for t in range(M):
	theta = t*2*pi/M
	circle[t] = np.array([np.cos(theta),np.sin(theta)])
#plt.plot(circle[:,0],circle[:,1],'r')
#plt.show() # Plotting the GLL-points and the circle

# Plotting the "Areals"

dist = (xis[1:]-xis[:N-1])/2
xis2 = xis[:N-1] + dist
one2 = np.ones(N-1)
D = np.zeros((N-1,N-1))

for i in range(N-1):

	plt.plot(xis2[i]*one2,xis2[:],'b')
	plt.plot(xis2[:],xis2[i]*one2,'b')

for i in range(N):
	for j in range(N):
		if (xis[i]**2+xis[j]**2<1):
			plt.plot(xis[i],xis[j],'r*')
			D[i,j] = 1


plt.plot(circle[:,0],circle[:,1],'r')
plt.show() # Plotting the GLL-points and the circle
areal_est = sum(sum(np.multiply(D,np.outer(dist,dist))))
print ('ESTIMATION DONE FOR {} GLL-POINTS'.format(N))
print ('estimated areal of circle:{}'.format(areal_est))
print ('actual areal of circle:{}'.format(pi/4))
print('coefficient:{}'.format(areal_est*4/pi))
