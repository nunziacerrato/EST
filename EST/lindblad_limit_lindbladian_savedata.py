''' In this program the limit Lindbladian is considered in order to save the following values:
    - Limit Lindbladian eigenvalues;
    - Entanglement Survival Time associated to each limit Lindbladian (this time is called t_ent);
    - The eigenvlues of the map Phi(t_ent).
'''
import time
import numpy as np
import pandas as pd
import tenpy
from lindblad import *
from lindblad_limit_lindbladian import *

# Parameters
N = 2
iterations = 20000
alpha, gamma = 1, 1
t_max = 40
t_step = 0.001
k_value = 'lim'

fixed_list = ['a']
common_path = "C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data"

# Compute the values of interest considering as Lindbladian the limit Lindbladian
t_in = time.time()
for fixed in fixed_list:
    dict_LindEigvals_tEnt = {'fixed' : fixed, 'Lind_eigvals' : [],\
                                        'Phi_tEnt_eigvals' : [], 't_ent': []}
    t_ent_list = []
    
    # For fixed k a given number (iterations) of Lindbladian matrices are constructed. For each
    # matrix, cycling on t with a fixed t_step, the negativity of entanglement is computed.
    # To achieve some information about the dependency of t_ent from the eigenvalues of the map
    # phi(t_ent), the eigenvalues of that superoperator at that time are also saved.
    for iteration in range(iterations):
        RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None))
        RM_H = tenpy.linalg.random_matrix.GUE((N,N))

        Lind_matr = Limit_Lindbladian_matrix_interact_pict(N,RM_D,RM_H,gamma)

        lind_eigvals = np.linalg.eigvals(Lind_matr)
        dict_LindEigvals_tEnt['Lind_eigvals'].append(list(lind_eigvals))

        t = 0
        while (t < t_max):
            n_ent = negat_ent(N,Lind_matr,t)
            if np.abs(np.real(n_ent)) < 1e-15:
                print(f'iter = {iteration}, t_ent = {t}')
                phit = phi_t(N,Lind_matr,t)
                phit_eigvals = np.linalg.eigvals(phit)
                dict_LindEigvals_tEnt['Phi_tEnt_eigvals'].append(list(phit_eigvals))
                dict_LindEigvals_tEnt['t_ent'].append(t)
                t_ent_list.append(t)
                break

            t += t_step
    
    df_tEnt = pd.DataFrame(t_ent_list, columns = ['t_ent'])

    # Costruct a DataFrame from the dictionary where all the values have been saved.
    # df_LindEigvals_tEnt = pd.DataFrame.from_dict(dict_LindEigvals_tEnt)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # !!!! To prevent ValueError: The array must be all the same lenght !!!!
    try:
        df_LindEigvals_tEnt = pd.concat([pd.Series(dict_LindEigvals_tEnt[key], name = f'{key}') \
                                            for key in dict_LindEigvals_tEnt.keys()], axis = 1)
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        

        # Save the DataFrame as an Excel file.
        # It will be possible to convert this file into a new DataFrame to analize data.
        writer = pd.ExcelWriter(f'{common_path}\\{fixed}_fixed_comparable_data\\Lindblad_eigvals_'\
            f't_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
        df_LindEigvals_tEnt.to_excel(writer, 'DataFrame')
        writer.save()
        print(f'FILE CORRECTLY SAVED FOR: k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}')
    except Exception as e:
        print(e)
        print(f'EXCEPTION FOUND FOR: k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}')
    finally:
        writer_emergency = pd.ExcelWriter(f'{common_path}\\{fixed}_fixed_comparable_data\\Only_t_ent_Lindblad_eigvals_'\
            f't_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
        df_tEnt.to_excel(writer_emergency, 'DataFrame')
        writer_emergency.save()
    
elapsed_time = time.time() - t_in
print(f'elapsed_time = {elapsed_time} s')