#! /bin/bash

import numpy as np

def gordon_hall_grid(gamma1, gamma2, gamma3, gamma4, xis, etas):
    """This quite extensive function generates the mesh for a deformed geometry,
    given the boundary mappings (gammas).
    INPUT:
        gammas: Functions or function handles that describe the boundary.
        xis: GLL-points on the reference square in the x-direction.
        etas: GLL-points in the y-direction.
    OUTPUT: Mesh matrices X and Y.
    NOTE: the numbering of boundary sides are important. Be careful to use the following
    convention:
        gamma1: x(1,etas),
        gamma2: x(xis,1),
        gamma3: x(-1,etas),
        gamma4: x(xis, -1).
    """
    Nx = len(xis)
    Ny = len(etas)

    #Interpolants:
    phi0 = lambda r: (1-r)/2.
    phi1 = lambda r: (1+r)/2.


    #constructing A-part:
    xA = np.zeros( (Ny,Nx) )
    yA = np.zeros( (Ny,Nx) )

    xA[:,Nx-1], yA[:,Nx-1] = gamma1(etas)
    xA[:,0], yA[:,0] = gamma3(etas)

    for i in range(Ny):
        xA[i,:] = (1-xis)/2.*xA[i,0] + (1+xis)/2.*xA[i,-1]
        yA[i,:] = (1-xis)/2.*yA[i,0] + (1+xis)/2.*yA[i,-1]

    #constructing the B-part:
    xB = np.zeros( (Ny,Nx) )
    yB = np.zeros( (Ny,Nx) )

    xB[0,:], yB[0,:] = gamma4(xis)
    xB[-1,:], yB[-1,:] = gamma2(xis)

    for i in range(Nx):
        xB[:,i] = (1-etas)/2.*xB[0,i] + (1+etas)/2.*xB[-1,i]
        yB[:,i] = (1-etas)/2.*yB[0,i] + (1+etas)/2.*yB[-1,i]

    #Constructing the C-part:
    xC = np.zeros( (Ny,Nx) )
    yC = np.zeros( (Ny,Nx) )

    for m in range(Ny):
        for n in range(Nx):
            xC[m,n] = (phi0(etas[m])*(phi0(xis[n])*xB[0,0] + phi1(xis[n])*xB[0,-1]) +
                    phi1(etas[m])*(phi0(xis[n])*xB[-1,0] + phi1(xis[n])*xB[-1,-1]) )
            yC[m,n] = (phi0(etas[m])*(phi0(xis[n])*yB[0,0] + phi1(xis[n])*yB[0,-1]) +
                    phi1(etas[m])*(phi0(xis[n])*yB[-1,0] + phi1(xis[n])*yB[-1,-1]) )

    #Total assembly and plot:
    X = xA + xB - xC
    Y = yA + yB - yC
    return X, Y
