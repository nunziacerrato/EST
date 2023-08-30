''' In this program the histograms of the entanglement surivival times are computed.
    Then, the rescaled histograms are plotted.
    Fixed 'a' or 'g' and different values of k are considered.
    '''

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from lindblad import *

# Parameters:
N = 3
iterations = 20000
alpha, gamma = 1, 1
median_ref = 1

fixed_list = ['g']
fixed_val = 1
k_list_alpha = [0.1,0.3,0.5,0.7,1,1.5,2]#,3,5,7,10]
k_list_gamma = [5,10,50,100,250,500,1000]
common_path = "C:\\Users\\cerra\\Documents\\GitHub\\EST_data_Ng2"

############################################ Histograms ############################################

# Histograms for fixed 'a' or 'g'
if True:
    # k_list_gamma = [5,10,50,100,250,500,1000]
    for fixed in fixed_list:
        if fixed == 'a':
            greek_fixed = r'$\alpha$'
            k_list = k_list_alpha
        if fixed == 'g':
            greek_fixed = r'$\gamma$'
            k_list = k_list_gamma
        
        fig_hist, ax_hist = plt.subplots(figsize=(15,10))
        # Slide title
        # ax_hist.set_title(fr'$\tau_{{ent}}$ distributions, '\
        #                   fr'N = {N}, {greek_fixed} = {fixed_val} fixed, {iterations} iterations',\
        #                   fontsize=35)
        # Paper title
        ax_hist.set_title(r'$\bf{P}$$_{ppt}^{(k)}(x)$, '\
                          fr'N = {N}, {iterations} iterations',\
                          fontsize=35)

        for k_value in k_list:
            if fixed == 'a':
                if k_value == 0:
                    print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                    break
                else:
                    alpha = fixed_val
                    gamma = alpha/k_value
            if fixed == 'g':
                gamma = fixed_val
                alpha = gamma*k_value
            
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}={fixed_val}_fixed\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present
            t_ent_values = df_fromexcel['t_ent'].dropna().values
            # Ccompute mean, standard deviation and median of t_ent values
            mean_t_ent = np.round(np.mean(t_ent_values),2)
            std_t_ent = np.round(np.std(t_ent_values),2)
            median_t_ent = np.round(np.median(t_ent_values),2)
            min_t_ent = np.round(min(t_ent_values),2)

            
            # Compute the histogram of t_ent values
            ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                    density = True, label = fr'k = {k_value}')
                    # , $\alpha$ = {np.round(alpha,2)}, '\ fr'$\gamma$ = {np.round(gamma,2)}')
            ax_hist.tick_params(labelsize=20)
            # Slide axes labels
            # ax_hist.set_xlabel(fr'$\tau_{{ent}}$', fontsize = 30)
            # ax_hist.set_ylabel(fr'Density', fontsize = 30)
            
            # Paper axes labels
            ax_hist.set_xlabel(fr'$x$', fontsize = 30)
            ax_hist.legend(fontsize = 22)
            fig_hist.savefig(f'{common_path}\\{fixed}={fixed_val}_fixed_plot\\Limit_Histogram_{fixed_list[0]}_fixed_N={N}', bbox_inches='tight')


    # plt.show()


############################################## Check ##############################################

# Prove that t_ent values for fixed a=a' (t') can be obtained from t_ent values for fixed g=g' (t)
# using the relation t = t'*((k g')/a'). In particular, prove that the distributions obtained
# numerically for fixed values of a coincide with the ones obtained by rescaling t_ent with the
# aforementioned relation.
if False:
    k_list = [0.1,0.3,0.5,0.7,1,1.5,2]
    fixed_list = ['a','g']
    greek_fixed_list = [r'$\alpha$',r'$\gamma$']
    fixed_val = 1
    fig_hist, ax_hist = plt.subplots(figsize=(15,10))
    
    for fixed, greek_fixed in zip(fixed_list,greek_fixed_list):
        ax_hist.set_title(fr'Comparison between numerical and rescaled $\tau_{{ent}}$ distributions,'\
                          fr' N = {N}, {greek_fixed} = {fixed_val} fixed, {iterations} iterations')
        for k_value in k_list:
            if fixed == 'a':
                alpha = fixed_val
                gamma = alpha/k_value
                # Create a DataFrame from the excel file where all the data are saved
                df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
                f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
                
                # Access to the column 't_ent' and drop NaN values if present
                t_ent_values = df_fromexcel['t_ent'].dropna().values

                # Compute the histogram of t_ent values
                ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                        density = True, label = fr'Numerical - k = {k_value}, '\
                        fr'$\alpha$ = {np.round(alpha,2)}, '\
                        fr'$\gamma$ = {np.round(gamma,2)}')
                ax_hist.set_xlabel(\
                                fr'$\tau_{{ent}}$')
                ax_hist.legend()
            
            if fixed == 'g':
                gamma = fixed_val
                alpha = gamma*k_value
                # Create a DataFrame from the excel file where all the data are saved
                df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
                f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
                
                # Access to the column 't_ent' and drop NaN values if present
                t_ent_values = (df_fromexcel['t_ent'].dropna().values)*k_value
                
                # Compute the histogram of t_ent values
                ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                        density = True, label = fr'Rescaled  - k = {k_value}, '\
                        fr'$\alpha$ = {np.round(alpha,2)}, '\
                        fr'$\gamma$ = {np.round(gamma,2)}')
                ax_hist.set_xlabel(\
                                fr'$\tau_{{ent}}$')
                ax_hist.legend()
    plt.show()

# Histograms - only Dissipator - for fixed 'g' and g = 1, g = 2
if False:
    fixed = 'g'
    greek_fixed = r'$\gamma$'
    k_value = 0
    alpha = 0
    fixed_val_list = [1]

    fig_hist, ax_hist = plt.subplots(figsize=(15,10))
    ax_hist.set_title(fr'$\tau_{{ent}}$ distributions, '\
                          fr'N = {N}, $\alpha$ = 0, {greek_fixed} fixed, {iterations} iterations')   
    
    for fixed_val in fixed_val_list:
        gamma = fixed_val 
        # Create a DataFrame from the excel file where all the data are saved
        df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
        f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
        
        # Access to the column 't_ent' and drop NaN values if present
        t_ent_values = df_fromexcel['t_ent'].dropna().values
        # Ccompute mean, standard deviation and median of t_ent values
        mean_t_ent = np.round(np.mean(t_ent_values),2)
        std_t_ent = np.round(np.std(t_ent_values),2)
        median_t_ent = np.round(np.median(t_ent_values),2)
        min_t_ent = np.round(min(t_ent_values),2)

        
        # Compute the histogram of t_ent values
        ax_hist.hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
                fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}, '\
                fr'min($\tau_{{ent}}$) = {min_t_ent}')
        ax_hist.set_xlabel(\
                        fr'$\tau_{{ent}}$')
        ax_hist.legend()
    # plt.show()


############################################ Rescaling ############################################

# Rescaled histograms for fixed 'a' or 'g': t --> t/coeff with coeff = Me(t)/Me_ref;
# The histograms are rescaled so that their median is equal to one.
if False:
    k_list_alpha = [0.001,0.1,0.3,0.5,0.7,1,1.5,2]
    k_list_gamma = [0,0.1,0.3,0.5,0.7,1,1.5,5,10,50,100,250,500,1000]
    for fixed in fixed_list:
        if fixed == 'a':
            greek_fixed = r'$\alpha$'
            k_list = k_list_alpha
        if fixed == 'g':
            greek_fixed = r'$\gamma$'
            k_list = k_list_gamma
        
        fig_hist_resc, ax_hist_resc = plt.subplots(figsize=(15,10))
        ax_hist_resc.set_title(fr'Rescaled $\tau_{{ent}}$ distributions, N = {N}, '\
                    fr'{greek_fixed} fixed, {iterations} iterations', fontsize = 20)
        
        for k_value in k_list:
            if fixed == 'a':
                if k_value == 0:
                    print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                    break
                else:
                    alpha = fixed_val
                    gamma = alpha/k_value
            if fixed == 'g':
                gamma = fixed_val
                alpha = gamma*k_value
            
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
                f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present:
            t_ent_values = df_fromexcel['t_ent'].dropna().values
            # Ccompute mean, standard deviation and median of t_ent values
            mean_t_ent = np.mean(t_ent_values)
            devst_t_ent = np.std(t_ent_values)
            median_t_ent =np.median(t_ent_values)

            # Known the median, compute the rescaling coefficient to rescale all t_ent values in
            # order to make them have a median equal to 1.
            rescaling_coeff_median = median_t_ent/median_ref
            t_ent_median_rescaled = t_ent_values/rescaling_coeff_median

            # Compute mean, standard deviation and median of rescaled t_ent values
            mean_t_ent_resc = np.round(np.mean(t_ent_median_rescaled),2)
            std_t_ent_resc = np.round(np.std(t_ent_median_rescaled),2)
            median_t_ent_resc = np.round(np.median(t_ent_median_rescaled),2)
            
            # Compute the histogram of these rescaled t_ent values
            ax_hist_resc.hist(t_ent_median_rescaled, bins = 'auto', histtype='step', fill = False,\
                        density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                        fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                        fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent_resc}, $\sigma(\tau_{{ent}})$'\
                        fr'= {std_t_ent_resc}, Me($\tau_{{ent}}$) = {median_t_ent_resc}')
            ax_hist_resc.set_xlabel(fr'$\tau_{{ent}}/Me(\tau_{{ent}})$', fontsize = 18)
            ax_hist_resc.set_ylabel(fr'Density', fontsize = 18)
            ax_hist_resc.legend(fontsize = 12)
    # plt.show()

# Normalized histograms for fixed 'a' or 'g': t --> (t - t_min)/(t_max - t_min);
# The histograms are rescaled so that they all have t_min = 0 and t_max = 1.
if False:
    for fixed in fixed_list:
        if fixed == 'a':
            greek_fixed = r'$\alpha$'
            k_list = k_list_alpha
        if fixed == 'g':
            greek_fixed = r'$\gamma$'
            k_list = k_list_gamma
        
        fig_hist_resc, ax_hist_resc = plt.subplots(figsize=(15,10))
        fig_hist_resc.suptitle(fr'Normalized $\tau_{{ent}}$ distributions, N = {N}, '\
                    fr'{greek_fixed} fixed, {iterations} iterations')
        
        for k_value in k_list:
            if fixed == 'a':
                if k_value == 0:
                    print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                    break
                else:
                    alpha = fixed_val
                    gamma = alpha/k_value
            if fixed == 'g':
                gamma = fixed_val
                alpha = gamma*k_value
            
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
                f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present:
            t_ent_values = df_fromexcel['t_ent'].dropna().values
            # Ccompute mean, standard deviation and median of t_ent values
            mean_t_ent = np.mean(t_ent_values)
            std_t_ent = np.std(t_ent_values)
            median_t_ent =np.median(t_ent_values)
            
            # Compute the minimun and the maximum value of t_ent
            t_ent_min = min(t_ent_values)
            t_ent_max = max(t_ent_values)
            
            # Rescale all t_ent values
            t_ent_median_rescaled = (t_ent_values - t_ent_min)/(t_ent_max - t_ent_min)

            # Compute mean, standard deviation and median of rescaled t_ent values
            mean_t_ent_resc = np.round(np.mean(t_ent_median_rescaled),2)
            std_t_ent_resc = np.round(np.std(t_ent_median_rescaled),2)
            median_t_ent_resc = np.round(np.median(t_ent_median_rescaled),2)
            
            # Compute the histogram of these rescaled t_ent values
            ax_hist_resc.hist(t_ent_median_rescaled, bins = 'auto', histtype='step', fill = False,\
                        density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                        fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                        fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent_resc}, $\sigma(\tau_{{ent}})$'\
                        fr'= {std_t_ent_resc}, Me($\tau_{{ent}}$) = {median_t_ent_resc}')
            ax_hist_resc.set_xlabel(fr'$\tau_{{ent}}/[Me(\tau_{{ent}})/Me_{{ref}}]$')
            ax_hist_resc.legend()
    # plt.show()

# Standardized histograms for fixed 'a' or 'g': t --> (t - <t>)/std(t).
if False:
    fixed_list = ['g']
    fixed_val = 1
    for fixed in fixed_list:
        if fixed == 'a':
            greek_fixed = r'$\alpha$'
            k_list = k_list_alpha
        if fixed == 'g':
            greek_fixed = r'$\gamma$'
            k_list = k_list_gamma
        
        fig_hist_stand, ax_hist_stand = plt.subplots(figsize=(15,10))
        fig_hist_stand.suptitle(fr'Standardized $\tau_{{ent}}$ distributions, '\
                    fr'$\overline{{ \tau}}_{{ent}}$ = 0,  $\sigma(\tau_{{ent}})$ = 1 N = {N}, '\
                    fr'{greek_fixed} fixed, {iterations} iterations')
        
        for k_value in k_list:
            if fixed == 'a':
                if k_value == 0:
                    print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                    break
                else:
                    alpha = fixed_val
                    gamma = alpha/k_value
            if fixed == 'g':
                gamma = fixed_val
                alpha = gamma*k_value
            
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present:
            t_ent_values = df_fromexcel['t_ent'].dropna().values
            # Ccompute mean, standard deviation and median of t_ent values
            mean_t_ent = np.mean(t_ent_values)
            std_t_ent = np.std(t_ent_values)
            median_t_ent =np.median(t_ent_values)

            # Compute the standardized t_ent values.
            t_ent_standardized = (t_ent_values - mean_t_ent)/std_t_ent

            # Compute mean, standard deviation and median of rescaled t_ent values
            mean_t_ent_stand = np.round(np.mean(t_ent_standardized),2)
            std_t_ent_stand = np.round(np.std(t_ent_standardized),2)
            median_t_ent_stand = np.round(np.median(t_ent_standardized),2)
            
            # Compute the histogram of these rescaled t_ent values
            ax_hist_stand.hist(t_ent_standardized, bins = 'auto', histtype='step', fill = False,\
                    density = True, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                    fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                    fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent_stand}, $\sigma(\tau_{{ent}})$'\
                    fr'= {std_t_ent_stand}, Me($\tau_{{ent}}$) = {median_t_ent_stand}')
            ax_hist_stand.set_xlabel(\
                            fr'$(\tau_{{ent}}-\overline{{ \tau}}_{{ent}})/\sigma(\tau_{{ent}})$')
            ax_hist_stand.legend()
    # plt.show()

# Compare histograms for the same fixed k
if False:

    k_list = [0.1,0.5,0.7,1.5]
    fig_hist, ax_hist = plt.subplots(figsize=(15,10), nrows = 2, ncols = 2, sharex = True)
    fig_hist.suptitle(fr'Compare $\tau_{{ent}}$ distributions, '\
                      fr'N = {N}, {iterations} iterations')
    
    for k_value in k_list:
        x_index = k_list.index(k_value)//2
        y_index = k_list.index(k_value)%2
        for fixed in ['a','g']:
            if fixed == 'a':
                color_ = 'orange'
                if k_value == 0:
                    print('k CANNOT BE EQUAL TO 0 IF ALPHA IS FIXED.')
                    break
                else:
                    alpha = fixed_val
                    gamma = alpha/k_value
            if fixed == 'g':
                color_ = 'blue'
                gamma = fixed_val
                alpha = gamma*k_value
            # Create a DataFrame from the excel file where all the data are saved
            df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')
            
            # Access to the column 't_ent' and drop NaN values if present
            t_ent_values = df_fromexcel['t_ent'].dropna().values
            # Ccompute mean, standard deviation and median of t_ent values
            mean_t_ent = np.round(np.mean(t_ent_values),2)
            std_t_ent = np.round(np.std(t_ent_values),2)
            median_t_ent = np.round(np.median(t_ent_values),2)

            
            # Compute the histogram of t_ent values
            ax_hist[x_index,y_index].hist(t_ent_values, bins = 'auto', histtype='step', fill = False,\
                    density = True, color = color_, label = fr'k = {k_value}, $\alpha$ = {np.round(alpha,2)}, '\
                    fr'$\gamma$ = {np.round(gamma,2)}'+ '\n' + \
                    fr'$\overline{{ \tau}}_{{ent}}$ = {mean_t_ent}, $\sigma(\tau_{{ent}})$'\
                    fr'= {std_t_ent}, Me($\tau_{{ent}}$) = {median_t_ent}')
            ax_hist[x_index,y_index].set_xlabel(fr'$\tau_{{ent}}$')
            ax_hist[x_index,y_index].legend()


# plt.show()