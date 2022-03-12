''' This program creates the limit Lindbladian in interaction picture as a superoperator acting on an 
    input matrix and then computes it matrix expression. This code should be intended as a library.
    The limit Lindbladian is defined as the expression of the Lindbladian where the alpha parameter
    tends to infinity.
'''
import math
import numpy as np
import pandas as pd
from lindblad import *

######################################### RHO CONST ###############################################
################################### ZERO ORDER CONTRIBUTION #######################################

def Limit_Lindbladian_interact_pict(N,RM_D,RM_H,matrix,gamma):
    ''' Function that creates the limit Lindbladian in interaction picture as a superoperator acting
        on the input matrix, starting from the Kossakowki matrix, constructed from a random matrix
        sampled from the Ginibre ensemble, and an Hamiltonian matrix sampled from a GUE.
        Here it is ensured that the trace of the Kossakowski matrix is equal to N.
        The limit Lindbladian is intended as the Lindbladian in which it is performed an integral
        average over alpha from 0 to an alpha_cutoff.
        Parameters: N : int
                        Dimension of the input matrix.
                    RM_D : ndarray
                           Random matrix sampled from the Ginibre ensemble.
                           This matrix can be obtained using the QuTip library in the following way:
                           RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None)).
                    RM_H : ndarray
                           Random hamiltonian matrix sampled from the GUE ensemble.
                           This matrix can be obtained using the tenpy library in the following way:
                           RM_H = tenpy.linalg.random_matrix.GUE((N,N)).
                    matrix : ndarray
                             The input matrix.
                    gamma : float
                            The coefficient that regulates the strenght of the dissipator.
        Returns:    out : ndarray
                          Output matrix after the action of a random limit lindbladian.
    '''
    K = N*RM_D

    # Diagonalize the normalized Wishart matrix -> Kossakowski matrix
    eigval_K, eigvect_K = np.linalg.eigh(K)

    # Diagonalize the Hamiltonian
    eigval_H, eigvect_H = np.linalg.eigh(RM_H)
    eigvect_H_dagg = (np.conjugate(eigvect_H)).T
    eigvect_H_inv = np.linalg.inv(eigvect_H)

    # Build Lindblad operators as an array of three indices: N*2 - 1 operators of dimension (N x N)
    F = F_matr_base_hs(N)

    L = np.zeros((N**2 -1,N,N), dtype=complex)
    # L_dagg = np.zeros((N**2 -1,N,N), dtype=complex)
    L_prime = np.zeros((N**2 -1,N,N), dtype=complex)

    I1_1 = np.zeros((N**2 -1,N,N), dtype=complex)
    I1_2 = np.zeros((N**2 -1,N,N), dtype=complex)
    I2 = np.zeros((N**2 -1,N,N), dtype=complex)
    I3 = np.zeros((N**2 -1,N,N), dtype=complex)

    Lind_lim = np.zeros((N,N), dtype=complex)

    # Compute the terms that made up the limit Lindbladian
    for k in range(N**2 -1):
        l = np.zeros((N,N), dtype=complex)
        for m in range(N**2 -1):
            l = l + eigvect_K[m,k]*F[m+1]  # You have to exclude the first element of F, Id(N).
        l = l*np.sqrt(eigval_K[k])
        L[k] = l
        # L_dagg[k] = (np.conjugate(l)).T

        L_prime[k] = eigvect_H_inv@(L[k])@eigvect_H

        for i in range(N):
            ii_outer = np.outer(eigvect_H[:,i],eigvect_H_dagg[i,:])
            Z_ii = eigvect_H@ii_outer@L_prime[k]@ii_outer@eigvect_H_inv
            for j in range(N):
                jj_outer = np.outer(eigvect_H[:,j],eigvect_H_dagg[j,:])
                Z_jj_dagg = eigvect_H@jj_outer@((np.conjugate(L_prime[k])).T)@jj_outer@eigvect_H_inv
                I1_1[k] = I1_1[k] + Z_ii@matrix@Z_jj_dagg

                Z_ji = eigvect_H@jj_outer@L_prime[k]@ii_outer@eigvect_H_inv
                Z_ji_dagg = np.conjugate(Z_ji).T

                I2[k] = I2[k] + Z_ji_dagg@Z_ji@matrix
                I3[k] = I3[k] + matrix@Z_ji_dagg@Z_ji

                if j!=i:
                    Z_ij = eigvect_H@ii_outer@L_prime[k]@jj_outer@eigvect_H_inv
                    Z_ij_dagg = np.conjugate(Z_ij).T
                    I1_2[k] = I1_2[k] + Z_ij@matrix@Z_ij_dagg
        
        Lind_lim = Lind_lim + gamma*(I1_1[k] + I1_2[k] - 0.5*I2[k] - 0.5*I3[k])
    
    return Lind_lim

def Limit_Lindbladian_matrix_interact_pict(N,RM_D,RM_H,gamma):
    ''' Function that calculates the matrix associated with Lindblad’s superoperator written with
        respect to the Hilbert-Schmidt matrix base. Called F[m] these matrices, for m = 1,...,N**2,
        the elements of the Lindbladian matrix are: L[m,m]=Tr(F[m]L(F[n])).
        Parameters: N : int
                        Dimension of the input matrix.
                    RM_D : ndarray
                           Random matrix sampled from the Ginibre ensemble and used to construct
                           the dissipator.
                           This matrix can be obtained using the QuTip library in the following way:
                           RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None)).
                    RM_H : ndarray
                           Random hamiltonian matrix sampled from the GUE ensemble.
                           This matrix can be obtained using the tenpy library in the following way:
                           RM_H = tenpy.linalg.random_matrix.GUE((N,N)).
                    gamma : float
                            The coefficient that regulates the strenght of the dissipator.
        Returns:    out : ndarray
                    Lindbladian matrix of dimension (N**2 x N**2), written in the Hilbert-Schmidt
                    matrix base.
    '''
    FF = F_matr_base_hs(N)
    lindbladian_matr = np.zeros((N**2,N**2), dtype=complex)

    for m in range(N**2):
        for n in range(N**2):
            A = FF[m]@Limit_Lindbladian_interact_pict(N,RM_D,RM_H,FF[n],gamma)
            lindbladian_matr[m,n] = np.trace(A)

    return lindbladian_matr
