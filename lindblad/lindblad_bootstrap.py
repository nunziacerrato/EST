''' In this program, for each value of k and for fixed value of alpha or gamma, new samples are
    created using the bootstrap method. The number of new saples is fixed at 100 and their dimension
    coincides with the number of points that one has in the inital sample.
    For each new sample the following quantities are saved: min(t_ent), Me(t_ent), mean(t_ent).
    So one has 100 of these quantities and it is possible to compute their mean value and their
    standard deviation. These data are finally saved in an Excel file.
'''
import numpy as np
import pandas as pd
import sklearn

# Parameters:
N = 2
iterations = 20000
alpha_fixed, gamma_fixed = 1, 1

bootstrap_sample_number = 100
bootstrap_sample_lenght = iterations

fixed_list = ['g']
fixed_val = 1
k_list_alpha = [0.1,0.3,0.5,0.7,1,1.5,2,3,5,7,10]
k_list_gamma = [0,0.1,0.3,0.5,0.7,1.5,5,10,50,100,250,500,1000]
k_list_gamma = [0,0.1,0.3,0.5,0.7,1,1.5,2,3,4,5,7,10,15,20,25,30,35,40,45,50,100,250,500,1000]

common_path = "C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data"

for fixed in fixed_list:
    if fixed == 'a':
        greek_fixed = r'$\alpha$'
        alpha = fixed_val
        k_list = k_list_alpha
    if fixed == 'g':
        greek_fixed = r'$\gamma$'
        k_list = k_list_gamma
        gamma = fixed_val

    bootstrap_mean_tent_estim = {'k' : [], 'alpha' : [], 'gamma' : [], '<min(t_ent)>' : [], 'min_std' : [],\
                          '<Me(t_ent)>' : [], 'Me_std' : [], '<mean(t_ent)>' : [], 'mean_std' : []}
    
    for k_value in k_list:
        if fixed == 'a':
            if k_value == 0:
                print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                break
            else:
                gamma = alpha/k_value
        if fixed == 'g':
            alpha = gamma*k_value
        
        # Create the dictionary where the values will be saved
        bootstrap_tent_estim = {'min(t_ent)' : [], 'Me(t_ent)' : [], 'mean(t_ent)' : []}


        # Create a DataFrame from the excel file where all the data are saved
        df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\t_step_more_precise_for_k_in_curvature\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
        
        
        # Access to the column 't_ent' and drop nan values:
        t_ent_values = df_fromexcel['t_ent'].dropna().values


        # Create new samples using the bootstrap method (see docs to set random_state as seed).
        # From the inital sample, create 100 new samples of the same dimension of the initial one
        # and save each sample in a dictionary, at the corresponding key.
        # For each sample extract the following values: min(t_ent), Me(t_ent), mean(t_ent), so there
        # will be 100 of these values. Compute the mean of these values and add them in a new dict.
        # Repeat this process 10 times, so there will be 10 of these means.
        # Finally compute the mean of means of these values and the stdev.

        # Generate 100 new samples (100 = bootstrap_sample_number)
        for sample_number in range(0,bootstrap_sample_number):
            sample = np.array(sklearn.utils.resample(t_ent_values, replace=True, \
                                                        n_samples=bootstrap_sample_lenght))

            # Save the min(t_ent), Me(t_ent), mean(t_ent) of the new sample
            bootstrap_tent_estim['min(t_ent)'].append(min(sample))
            bootstrap_tent_estim['Me(t_ent)'].append(np.median(sample))
            bootstrap_tent_estim['mean(t_ent)'].append(np.mean(sample))
        
        print(f'k = {k_value}, \
                bootstrap_tent_estim - must have {bootstrap_sample_number} elements')
        df_bootstrap_tent_estim = pd.DataFrame.from_dict(bootstrap_tent_estim)
        print(df_bootstrap_tent_estim)

        # Compute the mean values of the previous quantities and the standard deviation:
        # <min(t_ent)>, <Me(t_ent)>, <mean(t_ent)> and the standard deviation
        bootstrap_mean_tent_estim['k'].append(k_value)
        bootstrap_mean_tent_estim['alpha'].append(alpha)
        bootstrap_mean_tent_estim['gamma'].append(gamma)

        bootstrap_mean_tent_estim['<min(t_ent)>'].append(np.array( \
                                                    bootstrap_tent_estim['min(t_ent)']).mean())
        bootstrap_mean_tent_estim['min_std'].append(np.array( \
                                                    bootstrap_tent_estim['min(t_ent)']).std())
        bootstrap_mean_tent_estim['<Me(t_ent)>'].append(np.array( \
                                                    bootstrap_tent_estim['Me(t_ent)']).mean())
        bootstrap_mean_tent_estim['Me_std'].append(np.array( \
                                                    bootstrap_tent_estim['Me(t_ent)']).std())
        bootstrap_mean_tent_estim['<mean(t_ent)>'].append(np.array( \
                                                    bootstrap_tent_estim['mean(t_ent)']).mean())
        bootstrap_mean_tent_estim['mean_std'].append(np.array( \
                                                    bootstrap_tent_estim['mean(t_ent)']).std())
    
    print(f'df_mean_tent_estim - must have {len(k_list)} elements')
    df_mean_tent_estim = pd.DataFrame.from_dict(bootstrap_mean_tent_estim)
    print(df_mean_tent_estim)

    writer = pd.ExcelWriter(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\t_step_more_precise_for_k_in_curvature\\'\
        f'\\Lind_bootstrap_t_ent_{fixed}_fixed_{iterations}_iterations_N_{N}_more_precise_data_from_1_to_15.xlsx')
    df_mean_tent_estim.to_excel(writer, 'DataFrame')
    writer.save()
    print(f'FILE CORRECTLY SAVED FOR: k = {k_value}, a = {alpha}, g = {np.round(gamma,2)}')

