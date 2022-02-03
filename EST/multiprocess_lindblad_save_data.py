import time
import numpy as np
import pandas as pd
import multiprocessing
import tenpy
import qutip
from lindblad import *


# k_list_smallvalues = [0.1,0.3,0.5,0.7,1.5]
k_list_smallvalues = [2,5,10,50,100]
k_list_bigvalues = [1.5,2,2.5,3,5,7]

k_list_alpha = [0.1,0.3,0.5,0.7,1,1.5,2]
k_list_gamma = [0,0.1,0.3,0.5,0.7,1,1.5,5,10,50,100,250,500,1000]
k_list_new_terms_min5 = [2,2.5,3,3.5,4]
k_list_new_terms = [15,20,25,30,35,40,45]
k_list_more_precision_t_step = [1.5,2,3,4,5,7,10,15]

def lindblad_savedata(k_value):

    N = 2
    iterations = 20000
    fixed_val = 1
    alpha, gamma = 1, fixed_val
    t_max = 40
    t_step = 0.0001

    fixed_list = ['g']

    common_path = "C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data"

    # Iterations on k_values. Check if alpha or gamma are fixed and compute the other value.
    # The value k = 0 is taken into account when gamma is fixed, so it is possible to have alpha = 0,
    # gamma = 1. In the other case, when alpha is fixed, there could be an infinity computing the gamma
    # value when the value k = 0 in encountered, so in that case the algorithm passes.
    for fixed in fixed_list:
        if fixed == 'a':
            if k_value == 0:
                continue
            else:
                gamma = alpha/k_value
        if fixed == 'g':
            gamma = fixed_val
            if k_value == 1:
                continue
            alpha = gamma*k_value

        dict_LindEigvals_tEnt = {'k' : k_value, 'fixed' : fixed, 'Lind_eigvals' : [],\
                                            'Phi_tEnt_eigvals' : [], 't_ent': []}
        t_ent_list = []

        # For fixed k a given number (iterations) of Lindbladian matrices are constructed. For each
        # matrix, cycling on t with a fixed t_step, the negativity of entanglement is computed.
        # To achieve some informations about the dependency of t_ent from the eigenvalues of the map
        # phi(t_ent), the eigenvalues of that superoperator at that time are also saved.
        for iteration in range(iterations):
            RM_D = np.array(qutip.rand_dm_ginibre((N**2-1), rank=None))
            RM_H = tenpy.linalg.random_matrix.GUE((N,N))

            Lind_matr = Lindbladian_matrix(N,RM_D,RM_H,alpha,gamma)
            
            # lind_eigvals = np.linalg.eigvals(Lind_matr)
            # dict_LindEigvals_tEnt['Lind_eigvals'].append(list(lind_eigvals))
            
            # df_LindEigvals_tEnt = df_LindEigvals_tEnt.append({'Lind_eigvals': lind_eigvals},\
            #                                                 ignore_index=True)

            t = 0.78/gamma
            while (t < t_max):
                n_ent = negat_ent(N,Lind_matr,t)
                if np.abs(np.real(n_ent)) < 1e-13:
                    print(f'k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}, iter = {iteration}\
                            , t_ent = {t}')
                    # phit = phi_t(N,Lind_matr,t)
                    # phit_eigvals = np.linalg.eigvals(phit)
                    # dict_LindEigvals_tEnt['Phi_tEnt_eigvals'].append(list(phit_eigvals))
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
            # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            

            # Save the DataFrame as an Excel file.
            # It will be possible to convert this file into a new DataFrame to analize data.
            writer = pd.ExcelWriter(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\t_step_more_precise_for_k_in_curvature\\Lindblad_eigvals_'\
                f't_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            df_LindEigvals_tEnt.to_excel(writer, 'DataFrame')
            writer.save()
            print(f'FILE CORRECTLY SAVED FOR: k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}')
        except Exception as e:
            print(e)
            print(f'EXCEPTION FOUND FOR: k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}')
        finally:
            writer_emergency = pd.ExcelWriter(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\t_step_more_precise_for_k_in_curvature\\Only_t_ent_Lindblad_eigvals_'\
                f't_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            df_tEnt.to_excel(writer_emergency, 'DataFrame')
            writer_emergency.save()

if __name__== '__main__':
    t_in = time.time()
    processes = []
    k_list = k_list_more_precision_t_step
    for k_value in k_list:
        process = multiprocessing.Process(target=lindblad_savedata, 
                                          args=(k_value,))
        processes.append(process)
        process.start()
    for proc in processes:
        proc.join()
    print(f'elapsed time = {time.time() - t_in}')



    


