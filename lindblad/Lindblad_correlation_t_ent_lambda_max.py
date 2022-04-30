''' In this program we study the correlations between the inverse of the Lindbladian eigenvalue
    with the minimum real part, in modulo, excluding the null eigenvalue, and the entanglement
    survival time.
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.stats import pearsonr
import seaborn as sns
from lindblad import *

#####################################     DATA COLLECTION      #####################################
# Save the data on file: Lindbladian eigenvalues and entanglement survival time.
if False:
    N = 2
    alpha = 1
    gamma = 1
    k_list = [0.5,1,1.5,2,2.5,3,5,10]
    iterazioni = 2000
    path = 'C:\\Users\\prppg\\Desktop\\NUNZIA\\Dati\\Relaz_t_ent_autov_max'

    for k in k_list:
        index = 1
        alpha = k*gamma
        # gamma = alpha/k

        Lind_eigvals_txt = f"{path}\\Autovalori L per N = {N}, alpha = {alpha}, gamma = {gamma}, {iterazioni} iterazioni"+'.txt'
        f_Lind_eigvals = open(Lind_eigvals_txt, 'w')

        Lind_eigval_max_txt = f"{path}\\Autov max L per N = {N}, alpha = {alpha}, gamma = {gamma}, {iterazioni} iterazioni"+'.txt'
        f_Lind_eigval_max = open(Lind_eigval_max_txt, 'w')
        index = index + 1

        #################################################################
        
        t_ent_txt = f"C{path}\\t_ent per N = {N}, alpha = {alpha}, gamma = {gamma}, {iterazioni} iterazioni"+'.txt'
        f_t_ent = open(t_ent_txt, 'w')

        for iteraz in range(iterazioni):
            # Sampling the Lindbladian
            RM_D = np.array(qutip.rand_dm_ginibre((N**2 - 1), rank=None))
            RM_H = tenpy.linalg.random_matrix.GUE((N,N))

            Lind_matr = Lindbladian_matrix(N,RM_D,RM_H,alpha,gamma)

            # Writing the eigenvalues of L and the maximum eigenvalue of L in two different files
            eigvals_L = np.linalg.eigvals(Lind_matr)
            f_Lind_eigvals.write(f'{eigvals_L} \n')

            for item in range(N**2):
                if np.abs(np.real(eigvals_L[item])) < 1e-12:
                    eigvals_L_without0 = list(np.delete(eigvals_L,item))
            
            max_eigval_L = np.abs(np.real(np.max(eigvals_L_without0)))
            f_Lind_eigval_max.write(f'{max_eigval_L} \n')         
            
            # Writin in a file the time in which t_ent becomes null for that particular Lindbladian
            t = 0.
            t_max = 20.
            while(t<t_max):
                n_ent = negat_ent(N,Lind_matr,t)
                if np.abs(np.real(n_ent)) < 1e-15:
                    f_t_ent.write(f'{t} \n')
                    print(f'k = {k}, iteraz = {iteraz}, t = {t}, negat_ent = {n_ent}')
                    break
                t = t + 0.001


        f_Lind_eigvals.close()
        f_Lind_eigval_max.close()
        f_t_ent.close()


##########################################      PLOT      ##########################################
if True:
    N = 2
    k_list = [0.5,1,2,5]
    iterazioni = 2000

    alpha = 1
    gamma = 1
    path = 'C:\\Users\\cerra\\Desktop\\Thesis_data_plot\\Dati\\Relaz_t_ent_autov_max'

    for k in k_list:
        alpha = k*gamma
        # gamma = alpha/k

        # Writing in a list the maximum eigenvalue of L 
        Lind_eigval_max_txt = f"{path}\\Autov max L per N = {N}, alpha = {alpha}, gamma = {gamma}, {iterazioni} iterazioni.txt"
        f_Lind_eigval_max = open(Lind_eigval_max_txt, 'r')
        eigval_max = f_Lind_eigval_max.read()
        eigval_max_list = list(np.fromstring(eigval_max, sep = '\n'))

        # Computing the inverse of the maximum eigenvalue of L
        inverse_eigval_max_list = []
        for item in eigval_max_list:
            inverse_eigval_max_list.append(1/item)

        f_Lind_eigval_max.close()


        # Writing in a list the entanglement survival time, t_ent, corresponding to the Lindbladian\
        # of which I know the inverse of the maximum eigenvalue
        t_ent_txt = f"{path}\\t_ent per N = {N}, alpha = {alpha}, gamma = {gamma}, {iterazioni} iterazioni.txt"
        f_t_ent = open(t_ent_txt, 'r')
        t_ent = f_t_ent.read()
        t_ent_list = list(np.fromstring(t_ent, sep = '\n'))

        f_t_ent.close()
        
        correlation_pear = pearsonr(inverse_eigval_max_list, t_ent_list)
        correlation_spear = spearmanr(inverse_eigval_max_list, t_ent_list)
        print(f'correlation_pearson = {correlation_pear}')
        print(f'correlation_spearman = {correlation_spear}')


        p1 = sns.jointplot(x=inverse_eigval_max_list, y=t_ent_list, kind='reg')
        p1.fig.suptitle(fr'$\tau_{{ent}}$ as a function of 1/|$\lambda_{{max}}$|'\
            fr' - k={k}', fontsize=35)
        p1.ax_joint.legend([f'Spearman correlation = {np.round(correlation_spear[0],2)}; p-value = {correlation_spear[1]:.1E}'],\
            fontsize=30)

        p1.fig.subplots_adjust(top=0.93)

    
        plt.show()
