''' In this program the histograms are plotted for fixed 'g' and different values of k, together
    with the histogram obtained in correspondence of the limit Lindbladian, so for k tending to inf.
'''
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lindblad import *

# Parameters:
N = 2
iterations = 20000
alpha, gamma = 1, 1

fixed, greek_fixed = 'g', r'$\gamma$'
fixed_val = 1
k_list = [0,0.1,0.3,0.5,0.7,1.5,5,10,50,100]

common_path = "C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data"

# Histograms for fixed 'g' and different values of 'a'
if True:        
    fig_hist, ax_hist = plt.subplots(figsize=(15,10))
    ax_hist.set_title(fr'$\tau_{{ent}}$ distributions, '\
                        fr'N = {N}, {greek_fixed} = {fixed_val} fixed, {iterations} iterations',\
                        fontsize = 18)

    for k_value in k_list:
        gamma = fixed_val
        alpha = gamma*k_value
        
        # Create a DataFrame from the excel file where all the data are saved
        df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
        f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
        
        # Access to the column 't_ent' and drop NaN values if present
        t_ent_values = df_fromexcel['t_ent'].dropna().values
        # Ccompute mean, standard deviation and median of t_ent values
        mean_t_ent = np.round(np.mean(t_ent_values),3)
        std_t_ent = np.round(np.std(t_ent_values),3)
        median_t_ent = np.round(np.median(t_ent_values),3)
        min_t_ent = np.round(min(t_ent_values),3)
     
        # Compute the histogram of t_ent values
        ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
                fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}, '\
                fr'min($\tau_{{ent}}$) = {min_t_ent}')
        ax_hist.set_xlabel(\
                        fr'$\tau_{{ent}}$', fontsize = 16)
    
    df_fromexcel_lim = pd.read_excel(f'{common_path}\\a_fixed_comparable_data\\'\
                            f'Lindblad_eigvals_t_ent_a_fixed_k_lim_20000_iterations_N_2.xlsx')

    # Access to the column 't_ent' and drop NaN values if present
    t_ent_values_lim = df_fromexcel_lim['t_ent'].dropna().values
    mean_t_ent = np.round(np.mean(t_ent_values_lim),3)
    std_t_ent = np.round(np.std(t_ent_values_lim),3)
    median_t_ent = np.round(np.median(t_ent_values_lim),3)
    min_t_ent = np.round(min(t_ent_values_lim),3)
    ax_hist.hist(t_ent_values_lim, bins = 'auto', histtype='step', fill = False,\
                density = True, label = fr'Lind_lim $\alpha$ = inf, '\
                fr'$\gamma$ = 1'+ '\n' + \
                fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
                fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}, '\
                fr'min($\tau_{{ent}}$) = {min_t_ent}')
    ax_hist.legend()
    plt.show()
    # plt.savefig(f'C:\\Users\\cerra\\Desktop\\Problems_LL\\t_ent_distrib_with_lind_lim')

# Close up: histogram for the maximum value of 'a' considered and limit histogram
if True:
    k_value = 1000
    gamma = fixed_val

    fig_hist, ax_hist = plt.subplots(figsize=(15,10))
    ax_hist.set_title(fr'$\tau_{{ent}}$ distributions, '\
                      fr'N = {N}, {greek_fixed} = {fixed_val} fixed, {iterations} iterations',\
                      fontsize = 18)

    # Create a DataFrame from the excel file where all the data are saved
    df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
    f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
    
    # Access to the column 't_ent' and drop NaN values if present
    t_ent_values = df_fromexcel['t_ent'].dropna().values
    # Ccompute mean, standard deviation and median of t_ent values
    mean_t_ent = np.round(np.mean(t_ent_values),3)
    std_t_ent = np.round(np.std(t_ent_values),3)
    median_t_ent = np.round(np.median(t_ent_values),3)
    min_t_ent = np.round(min(t_ent_values),3)

    # Compute the histogram of t_ent values
    ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
            density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
            fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
            fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
            fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}, '\
            fr'min($\tau_{{ent}}$) = {min_t_ent}')
    ax_hist.set_xlabel(\
                    fr'$\tau_{{ent}}$', fontsize = 16)
        
    df_fromexcel_lim = pd.read_excel(f'{common_path}\\a_fixed_comparable_data\\'\
                                    f'Lindblad_eigvals_t_ent_a_fixed_k_lim_20000_iterations_N_2.xlsx')

    # Access to the column 't_ent' and drop NaN values if present
    t_ent_values_lim = df_fromexcel_lim['t_ent'].dropna().values
    mean_t_ent = np.round(np.mean(t_ent_values_lim),3)
    std_t_ent = np.round(np.std(t_ent_values_lim),3)
    median_t_ent = np.round(np.median(t_ent_values_lim),3)
    min_t_ent = np.round(min(t_ent_values_lim),3)
    ax_hist.hist(t_ent_values_lim, bins = 'auto', histtype='step', fill = False,\
                density = True, label = fr'Lind_lim $\alpha$ = inf, '\
                fr'$\gamma$ = 1'+ '\n' + \
                fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
                fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}, '\
                fr'min($\tau_{{ent}}$) = {min_t_ent}')
    ax_hist.legend()
    plt.show()
