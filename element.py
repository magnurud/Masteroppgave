from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import numpy as np
import matplotlib.pyplot as plt
import gh
import quadrature_nodes as qn
import structured_grids as sg 
import laplace_functions as lp 
#Number of GLL points:
N = 10
xis = qn.GLL_points(N)
weights = qn.GLL_weights(N, xis)
etas = xis

def plotjacobi(gamma1,gamma2,gamma3,gamma4):
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
    ##
    print "Assembling stiffness matrix."
    A = lp.assemble_local_stiffness_matrix(D, G_tot, N, weights)
        Jmean = np.mean(Jac)
        for i in range(N):
            #print np.max(Jac[N*i:N*(i+1)])-np.min(Jac[N*i:N*(i+1)])
            print np.max(Jac[N*i:N*(i+1)-1]-Jac[N*i+1:N*(i+1)])
            #print np.max(Jac[range(i,N*(N-2)+i,N)]-Jac[range(N+i,N*(N-1)+i,N)])
        #print np.linalg.cond(A)
        #print G_tot
        print np.max(Jac)-np.min(Jac)
        ##Plot all rows:
        for i in range(N):
                        plt.plot(X[i,:],Y[i,:],'b')
                        plt.plot(X[:,i],Y[:,i],'b')
        fig = plt.figure()

        ax = fig.gca(projection='3d')
        surf = ax.plot_surface(X, Y, (Jac-Jmean).reshape(N,N), rstride=1, cstride=1, cmap=cm.coolwarm,
                linewidth=0, antialiased=False)
        plt.show()


        # CHECKING THE BOUNDARIES
        x = np.zeros( (N,4) ) #Checking the boundaries
        y = np.zeros( (N,4) ) #Checking the boundaries
        x[:,0],y[:,0] = gamma1(xis)
        x[:,1],y[:,1] = gamma2(xis)
    x[:,2],y[:,2] = gamma3(xis)
    x[:,3],y[:,3] = gamma4(xis)
    plt.plot(x[:,0],y[:,0],'b')
    plt.plot(x[:,1],y[:,1],'r')
    plt.plot(x[:,2],y[:,2],'g')
    plt.plot(x[:,3],y[:,3],'c')
    plt.show()
    return
