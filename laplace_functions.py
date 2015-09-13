# This module contains function to be used to assemble the stiffness matrix corresponding to the Laplacian.

import numpy as np


def G_matrix(xD, xTD, yD, yTD, alpha, beta):
    """Function returning the G-matrix for the stiffness-quadrature
    on the reference element [-1,1]^2.
    IMPORTANT UPDATE: Now also returns the Jacobian.
    INPUT:
        xD: The local x-matrix multiplied with the lagrangian derivative matrix.
        xTD: x-matrix transposed multiplied with D.
        yD, yTD: Analogous to xD and xTD.
        alpha, beta: Integer indices determining which point should be evaluated.

    OUTPUT:
        G: a 2-by-2 matrix giving the Jacobian and gradient transform in a point.
        G = J*Delta_N(xi_alpha, xi_beta)
        J: Jacobian at that point.
    """
    x_xi = xD[beta, alpha]
    x_eta = xTD[alpha, beta]
    y_xi = yD[beta, alpha]
    y_eta = yTD[alpha, beta]

    #Determining the Jacobian:
    J = x_xi*y_eta - x_eta*y_xi

    return J, 1./J * np.array( [ [x_eta**2 + y_eta**2, -(y_eta*y_xi + x_eta*x_xi)],
            [-(y_eta*y_xi + x_eta*x_xi), x_xi**2+y_xi**2] ])


def assemble_total_G_matrix(xD, xTD, yD, yTD, nx, ny):
    """Function calculating all the 2-by-2 G-matrices for all points in the element.
    IMPORTANT UPDATE: Now also returns the Jacobian at all local points.
    INPUT:
        xD: Local x-coordinate matrix multiplied with lagrangian derivative matrix.
        xTD: The transpose of the local x-matrix and the D-matrix.
        yD, yTD: Analogously for the y-coordinates.
    OUTPUT:
        G: total G-matrix for all points in the element.
    """
    # Number of points in the element:
    num_points = nx*ny

    #Initialise the G-matrix:
    G = np.zeros( (num_points, 2, 2) )
    Jac = np.zeros( num_points )

    for I in range(num_points):
        i = I%nx
        j = I/nx
        Jac[I], G[I] = G_matrix(xD, xTD, yD, yTD, i, j)

    return Jac, G


def calculate_local_stiffness_element( I, J, D, G_tot, weights,N ):
    """ Function for calculating the I,J element of the LOCAL
    stiffness matrix. Later to be assembled with the larger global
    stiffness matrix.
    INPUT:
        I,J: Local indices.
        D: Lagrangian derivative matrix.
        G_tot: G-matrices for all points.
        weights: GLL_weights.

    OUTPUT:
        A_K[I,J]: Element of the local stiffness matrix.
    """

    # Get x- and y- indices:
    i = I%N
    j = I/N

    k = J%N
    l = J/N

    #Initialise value:
    A = 0.
    A += weights[i]*weights[l]*D[j,l]*D[k,i]*G_tot[l*N+i,1,0]
    A += weights[k]*weights[j]*D[l,j]*D[i,k]*G_tot[j*N+k,0,1]
    
    if (j==l):
        for alpha in range(N):
            A += weights[alpha]*weights[l]*D[k,alpha]*D[i,alpha]*G_tot[l*N+alpha, 0, 0]

    if (i == k):
        for beta in range(N):
            A += weights[k]*weights[beta]*D[j,beta]*D[l,beta]*G_tot[beta*N+k,1,1]

    return A


def assemble_local_stiffness_matrix(D, G_tot, N, weights):
    """ Function for assembling the local stiffness matrix for one element
    prescribed by X and Y matrices.
    INPUT:
        D: Lagrangian derivative matrix.
        G_tot: Total G-matrix for the element (3D array)
        N: Number of GLL-points in each direction.
        weights: GLL-weights.
    OUTPUT:
        A_K: Local stiffness matrix for element K.
    """

    # Initialise local stiffness matrix:
    num_points = N**2
    A_K = np.zeros( (num_points, num_points) )

    for I in range(num_points):
        for J in range(I+1):
            A_K[I,J] += calculate_local_stiffness_element(I, J, D, G_tot, weights,N)
            #if (J < I):
            #    A_K[J,I] = A_K[I,J]

    return A_K + A_K.T - np.diag(A_K.diagonal())

def assemble_loading_vector(X, Y, f, Jac, weights):
    """Function assembling the local loading vector F.
    INPUT:
        X, Y: Local coordinate matrices.
        f: Function or function handle for the loading function.
        Jac: Jacobian at all local points.
        weights: GLL-weights.
    OUTPUT:
        F: Local loading vector.
    """
    #Initialise loading vector.
    N = len(weights)
    num_points = N**2
    F = np.zeros( num_points )

    for I in range(num_points):
        i = I/N
        j = I%N
        F[I] += weights[i]*weights[j]*f(X[i,j],Y[i,j])*Jac[I]

    return F
    
    
    #THIS FUNCTION MAY BE BULLSHIT
def calculate_local_stiffness_element2( I, J, D_x,D_y, G_tot, weights_x,weights_y,Nx,Ny ):
    """ Function for calculating the I,J element of the LOCAL
    stiffness matrix. Later to be assembled with the larger global
    stiffness matrix.
    INPUT:
        I,J: Local indices.
        D: Lagrangian derivative matrix.
        G_tot: G-matrices for all points.
        weights: GLL_weights.

    OUTPUT:
        A_K[I,J]: Element of the local stiffness matrix.
    """

    # Get x- and y- indices:
    i = I%Nx
    j = I/Ny

    k = J%Nx
    l = J/Ny

    #Initialise value:
    A = 0.
    #THIS ESPECIALLY MIGHT BE FULL OF SHIT
    A += weights_x[i]*weights_y[l]*D_y[j,l]*D_x[k,i]*G_tot[l*Ny+i,1,0]
    A += weights_x[k]*weights_y[j]*D_y[l,j]*D_x[i,k]*G_tot[j*Ny+k,0,1]
    
    if (j==l):
        for alpha in range(Nx):
            A += weights_x[alpha]*weights_y[l]*D_x[k,alpha]*D_y[i,alpha]*G_tot[l*Ny+alpha, 0, 0]

    if (i == k):
        for beta in range(Ny):
            A += weights_x[k]*weights_y[beta]*D_x[j,beta]*D_y[l,beta]*G_tot[beta*Ny+k,1,1]

    return A

#THIS MIGHT BE TOTAL BOLLOCKS
def assemble_local_stiffness_matrix2(Dx,Dy, G_tot, Nx,Ny, weightsx,weightsy):

    # Initialise local stiffness matrix:
    num_points = Nx*Ny
    A_K = np.zeros( (num_points, num_points) )

    for I in range(num_points):
        for J in range(I+1):
            A_K[I,J] += calculate_local_stiffness_element2(I, J, Dx, Dy, G_tot, weightsx,weightsy,Nx,Ny)
            #if (J < I):
            #    A_K[J,I] = A_K[I,J]

    return A_K + A_K.T - np.diag(A_K.diagonal())


