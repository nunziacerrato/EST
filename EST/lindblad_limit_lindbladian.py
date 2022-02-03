import time
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tenpy
import qutip
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



if __name__ == '__main__':
    
    # Set parameters:
    N = 2
    gamma = 1
    t_max = 2.0
    t_step = 0.001
    common_path = 'C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data\\limit_lindbladian'
    
    for iter_ in range(20):
        print(f'iter = {iter_}')
        counts = 0
        # Plot of neg_ent(t) - Lindbladian for different values of alpha:
        if True:
            # Set random matrices
            RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None))
            RM_H = tenpy.linalg.random_matrix.GUE((N,N))

            # Plot of neg_ent(t)
            fig, ax = plt.subplots(figsize=(15,10))
            alpha_list = [0.001,0.01,0.1,1,3,5,10,100,1e3,1e4,1e5,1e6]#,1e7,1e8,1e9,1e10]
            t_ent_list = []
            for alpha in alpha_list:
                Lind_matr = Lindbladian_matrix(N,RM_D,RM_H,alpha,gamma)
                neg_ent_list = []
                t_list = []
                t = 0
                counts = 0
                while (t<t_max):
                    neg_ent = negat_ent(N,Lind_matr,t)
                    neg_ent_list.append(neg_ent)
                    t_list.append(t)
                    if np.abs(np.real(neg_ent)) < 1e-12 and alpha == alpha_list[-1] and counts == 0:
                        t_ent_alpha_max = t
                        print(f'alpha = {alpha:.1E}, t_ent = {t}')
                        ax.text(1.75, 0.23, fr'$\tau_{{ent}}$ for $\alpha=1e6$ = {np.round(t,3)}', fontsize=10)
                        counts = 1
                    t += t_step

                ax.plot(t_list, neg_ent_list, label = f'alpha = {alpha:.1E}')


            # Plot of neg_ent(t) with limit Lindbladian that does not depends on alpha
            if True:
                Lind_lim_matr = Limit_Lindbladian_matrix_interact_pict(N,RM_D,RM_H,gamma)
                neg_ent_lim_list = []
                t_lim_list = []
                # t_ent_lim = []
                t = 0
                counts = 0
                while (t<t_max):
                    neg_ent_lim = negat_ent(N,Lind_lim_matr,t)
                    neg_ent_lim_list.append(neg_ent_lim)
                    t_lim_list.append(t)
                    if np.abs(np.real(neg_ent_lim)) < 1e-12 and counts == 0:
                        t_ent_lim = t
                        print(f'alpha_lim, t_ent = {t}')
                        ax.text(1.75, 0.21, fr'$\tau_{{ent}}$ for Lind_lim = {np.round(t,3)}', fontsize=10)
                        counts = 1
                    t += t_step
                ax.plot(t_lim_list, neg_ent_lim_list, color = 'k', linestyle = 'dashed', label = f'Lind_limit')
        

                ax.set_title(fr'$neg_{{ent}} (t)$')
                ax.set_xlabel('t')
                ax.set_ylabel(fr'$neg_{{ent}}$')
                ax.legend()

                if t_ent_lim > t_ent_alpha_max:
                    print(f'error for iter = {iter_}')
                else:
                    print('good')

                plt.savefig(f'C:\\Users\\cerra\\Desktop\\Problems_LL\\new_LL_iter={iter_+1}')

    plt.show()



    # Percentili
    if False:
        alpha_list = [0,0.1,0.3,0.5,0.7,1,1.5,5,10,50,100,250,500]
        fixed = 'g'
        fixed_val = 1
        gamma = fixed_val
        iterations = 20000

        t_max = 15
        t_step = 0.001
        common_path = 'C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data'

        N = 2
        RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None))
        RM_H = tenpy.linalg.random_matrix.GUE((N,N))
        
        t_ent_list = []
        percentile_list = []
        for alpha in alpha_list:
            
            Lind_matr = Lindbladian_matrix(N,RM_D,RM_H,alpha,gamma)
            neg_ent_list = []
            t_list = []
            t = 0
            while (t<t_max):
                neg_ent = negat_ent(N,Lind_matr,t)
                neg_ent_list.append(neg_ent)
                t_list.append(t)
                if np.abs(np.real(neg_ent)) < 1e-12:
                    print(f'alpha = {alpha:.1E}, t_ent = {t}')
                    t_ent_list.append(t)
                    break
                t += t_step


            k_value = alpha
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present
            t_ent_values_in = list(df_fromexcel['t_ent'].dropna().values)
            t_ent_values = t_ent_values_in
            t_ent_values.sort()
            
            while len(t_ent_values)!=1:
                t_ent_values_central = t_ent_values[int(len(t_ent_values)/2)-1]
                if t <= t_ent_values_central:
                    t_ent_values = t_ent_values[:int(len(t_ent_values)/2)]
                else:
                    t_ent_values = t_ent_values[int(len(t_ent_values)/2):]
            
            index_ = t_ent_values_in.index(t_ent_values[0])
            percentile = index_/iterations
            percentile_list.append(percentile)
            print(f'alpha = {alpha}, percentile = {percentile}')
                
