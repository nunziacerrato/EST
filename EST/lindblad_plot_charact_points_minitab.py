''' In this program, the characteristics points of t_ent distributions are plotted as a function of
    alpha, having fixed gamma = 1. The points reported have been obtained using the Bootstrap method,
    and their error is 1 sigma. The fitted function is defined in the program and it is a ratio of
    polynomials, while the fit parameter that it takes in input have been obtained using MINITAB.
    For each graph, the curves obtained as upper and lower bound of the fit curve with respect
    to a 95% CI for the fit parameters are also reported.
'''
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# Parameters:
N = 2
iterations = 20000

common_path = 'C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data\\g_fixed_comparable_data\\g=1\\'\
                't_step_more_precise_for_k_in_curvature'

df_fromexcel = pd.read_excel(f'{common_path}\\'\
f'Lind_bootstrap_t_ent_g_fixed_{iterations}_iterations_N_{N}_more_precise_data_from_1_to_15.xlsx')

def fit_function(x,Theta1,Theta2,Theta3):
    return (Theta1*Theta3 + Theta2*(x**2))/(Theta3 + (x**2))


def min_t_ent_plot(x_values, k_max, plot_res = True):
    '''Plot the scatterplot of min_t_ent(alpha) until alpha = k_max with errors equals to 1 sigma,
    the fit function using as fit parameters the ones estimated with MINITAB, and the upper bound
    and lower bound fit functions obtained in correspondence of fit parameters that ensure a 95% CI.
    '''
    index_k_max = x_values.index(k_max)
    x_values = x_values[:index_k_max+1]
    min_t_ent = df_fromexcel['<min(t_ent)>'][:index_k_max+1]
    min_t_ent_error = df_fromexcel['min_std'][:index_k_max+1]

    plt.errorbar(x_values, min_t_ent, yerr = min_t_ent_error, fmt = 'o', markersize = 7,\
                label = fr'data  -  1$\sigma$ error')
    # fit parameters
    Theta1 = 0.837178
    Theta2 = 0.823919
    Theta3 = 0.697450
    # lower and upper bound of fit parameters - 95% CI
    Theta1_minus, Theta1_plus = 0.835772, 0.83864
    Theta2_minus, Theta2_plus = 0.823343, 0.82449
    Theta3_minus, Theta3_plus = 0.448919, 1.09220

    x_values_continue = np.linspace(x_values[0], x_values[-1], 500)

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1,Theta2,Theta3), label = \
        'fit function' + '\n'+ \
        fr'$\Theta_{{1}}$ = {Theta1}, $\Theta_{{2}}$ = {Theta2}, $\Theta_{{3}}$ = {Theta3}')
    
    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_minus,Theta2_minus,Theta3_minus),\
        color='green', linestyle='dashed', label = '95% CI')

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_plus,Theta2_plus,Theta3_plus),\
        color='green', linestyle='dashed')
    
    plt.xlabel(fr'$\alpha$', fontsize = 16)
    plt.ylabel(fr'$<min(\tau_{{ent}})>$', fontsize = 16)
    plt.title(fr'$<min(\tau_{{ent}})>(\alpha)$ with fit curve - 95% CI', fontsize = 18)
    plt.legend(fontsize = 14)

    if plot_res == True:
        
        fitted_values = fit_function(np.array(x_values),Theta1,Theta2,Theta3)
        residuals = min_t_ent - fitted_values

        fig_all, ax_all = plt.subplots(nrows=2, ncols=2)
        fig_all.suptitle(fr'Residual Plots for <min($\tau_{{ent}}$)>', fontsize = 18)

        # Plot Normal Probability Plot
        # fig_npp, ax_npp = plt.subplots()
        ax_npp = ax_all[0,0]
        res = stats.probplot(residuals, plot=ax_npp)

        # Plot residulas vs fitted values
        # fig_res, ax_res = plt.subplots()
        ax_res = ax_all[0,1]
        ax_res.scatter(fitted_values, residuals)
        ax_res.hlines(y=0, xmin=min(fitted_values), xmax=max(fitted_values), linestyle = 'dashed',\
                      color = 'black')
        ax_res.set_xlabel('Fitted Values')
        ax_res.set_ylabel('Residuals')
        ax_res.set_title('Versus Fit')
        # ax_res.set_xlim(min(fitted_values)-0.0001,max(fitted_values)+0.0001)
        
        # Plot histogram of residuals
        # fig_hist, ax_hist = plt.subplots()
        ax_hist = ax_all[1,0]
        ax_hist.hist(residuals, color='#0504aa', alpha=0.7, rwidth=0.9)
        ax_hist.set_xlabel('Residuals')
        ax_hist.set_ylabel('Counts')
        ax_hist.set_title('Histogram')

        # Plot residuals vs order
        # fig_order, ax_order = plt.subplots()
        ax_order = ax_all[1,1]
        order_list = np.arange(1,len(residuals)+1,1)
        ax_order.plot(order_list,residuals,'-o')
        ax_order.hlines(y=0, xmin=min(order_list), xmax=max(order_list), linestyle = 'dashed',\
                      color = 'black')
        ax_order.set_xlabel('Observation Order')
        ax_order.set_ylabel('Residuals')
        ax_order.set_title('Versus Order')


    plt.show()


def Me_t_ent_plot(x_values, k_max, plot_res = True):
    '''Plot the scatterplot of Me_t_ent(alpha) until alpha = k_max with errors equals to 1 sigma,
    the fit function using as fit parameters the ones estimated with MINITAB, and the upper bound
    and lower bound fit functions obtained in correspondence of fit parameters that ensure a 95% CI.
    '''
    index_k_max = x_values.index(k_max)
    x_values = x_values[:index_k_max+1]
    Me_t_ent = df_fromexcel['<Me(t_ent)>'][:index_k_max+1]
    Me_t_ent_error = df_fromexcel['Me_std'][:index_k_max+1]

    plt.errorbar(x_values, Me_t_ent, yerr = Me_t_ent_error, fmt = 'o', markersize = 7,\
                label = fr'data  -  1$\sigma$ error')
    # fit parameters
    Theta1 = 1.05137
    Theta2 = 0.85986
    Theta3 = 1.05223
    # lower and upper bound of fit parameters - 95% CI
    Theta1_minus, Theta1_plus = 1.04988, 1.05287
    Theta2_minus, Theta2_plus = 0.85923, 0.86049
    Theta3_minus, Theta3_plus = 1.01448, 1.09147

    x_values_continue = np.linspace(x_values[0], x_values[-1], 500)

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1,Theta2,Theta3), label = \
        'fit function' + '\n'+ \
        fr'$\Theta_{{1}}$ = {Theta1}, $\Theta_{{2}}$ = {Theta2}, $\Theta_{{3}}$ = {Theta3}')
    
    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_minus,Theta2_minus,Theta3_minus),\
        color='green', linestyle='dashed', label = '95% CI')

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_plus,Theta2_plus,Theta3_plus),\
        color='green', linestyle='dashed')
    
    plt.xlabel(fr'$\alpha$', fontsize = 16)
    plt.ylabel(fr'$<Me(\tau_{{ent}})>$', fontsize = 16)
    plt.title(fr'$<Me(\tau_{{ent}})>(\alpha)$ with fit curve - 95% CI', fontsize = 18)
    plt.legend(fontsize = 14)

    if plot_res == True:
        
        fitted_values = fit_function(np.array(x_values),Theta1,Theta2,Theta3)
        residuals = Me_t_ent - fitted_values

        fig_all, ax_all = plt.subplots(nrows=2, ncols=2)
        fig_all.suptitle(fr'Residual Plots for <Me($\tau_{{ent}}$)>', fontsize = 18)

        # Plot Normal Probability Plot
        # fig_npp, ax_npp = plt.subplots()
        ax_npp = ax_all[0,0]
        res = stats.probplot(residuals, plot=ax_npp)

        # Plot residulas vs fitted values
        # fig_res, ax_res = plt.subplots()
        ax_res = ax_all[0,1]
        ax_res.scatter(fitted_values, residuals)
        ax_res.hlines(y=0, xmin=min(fitted_values), xmax=max(fitted_values), linestyle = 'dashed',\
                      color = 'black')
        ax_res.set_xlabel('Fitted Values')
        ax_res.set_ylabel('Residuals')
        ax_res.set_title('Versus Fit')
        # ax_res.set_xlim(min(fitted_values)-0.0001,max(fitted_values)+0.0001)
        
        # Plot histogram of residuals
        # fig_hist, ax_hist = plt.subplots()
        ax_hist = ax_all[1,0]
        ax_hist.hist(residuals, color='#0504aa', alpha=0.7, rwidth=0.9)
        ax_hist.set_xlabel('Residuals')
        ax_hist.set_ylabel('Counts')
        ax_hist.set_title('Histogram')

        # Plot residuals vs order
        # fig_order, ax_order = plt.subplots()
        ax_order = ax_all[1,1]
        order_list = np.arange(1,len(residuals)+1,1)
        ax_order.plot(order_list,residuals,'-o')
        ax_order.hlines(y=0, xmin=min(order_list), xmax=max(order_list), linestyle = 'dashed',\
                      color = 'black')
        ax_order.set_xlabel('Observation Order')
        ax_order.set_ylabel('Residuals')
        ax_order.set_title('Versus Order')


    plt.show()


def mean_t_ent_plot(x_values, k_max, plot_res = True):
    '''Plot the scatterplot of mean_t_ent(alpha) until alpha = k_max with errors equals to 1 sigma,
    the fit function using as fit parameters the ones estimated with MINITAB, and the upper bound
    and lower bound fit functions obtained in correspondence of fit parameters that ensure a 95% CI.
    '''
    index_k_max = x_values.index(k_max)
    x_values = x_values[:index_k_max+1]
    mean_t_ent = df_fromexcel['<mean(t_ent)>'][:index_k_max+1]
    mean_t_ent_error = df_fromexcel['mean_std'][:index_k_max+1]

    plt.errorbar(x_values, mean_t_ent, yerr = mean_t_ent_error, fmt = 'o', markersize = 7,\
                label = fr'data  -  1$\sigma$ error')
    # fit parameters
    Theta1 = 1.08166
    Theta2 = 0.88504
    Theta3 = 0.97505
    # lower and upper bound of fit parameters - 95% CI
    Theta1_minus, Theta1_plus = 1.07948, 1.08384
    Theta2_minus, Theta2_plus = 0.88413, 0.88594
    Theta3_minus, Theta3_plus = 0.92584, 1.02705
    x_values_continue = np.linspace(x_values[0], x_values[-1], 500)

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1,Theta2,Theta3), label = \
        'fit function' + '\n'+ \
        fr'$\Theta_{{1}}$ = {Theta1}, $\Theta_{{2}}$ = {Theta2}, $\Theta_{{3}}$ = {Theta3}')
    
    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_minus,Theta2_minus,Theta3_minus),\
        color='green', linestyle='dashed', label = '95% CI')

    plt.plot(x_values_continue, fit_function(x_values_continue,Theta1_plus,Theta2_plus,Theta3_plus),\
        color='green', linestyle='dashed')
    
    plt.xlabel(fr'$\alpha$', fontsize = 16)
    plt.ylabel(fr'$<mean(\tau_{{ent}})>$', fontsize = 16)
    plt.title(fr'$<mean(\tau_{{ent}})>(\alpha)$ with fit curve - 95% CI',fontsize = 18)
    plt.legend(fontsize = 14)

    if plot_res == True:
        
        fitted_values = fit_function(np.array(x_values),Theta1,Theta2,Theta3)
        residuals = mean_t_ent - fitted_values

        fig_all, ax_all = plt.subplots(nrows=2, ncols=2)
        fig_all.suptitle(fr'Residual Plots for <mean($\tau_{{ent}}$)>',fontsize = 18)

        # Plot Normal Probability Plot
        # fig_npp, ax_npp = plt.subplots()
        ax_npp = ax_all[0,0]
        res = stats.probplot(residuals, plot=ax_npp)

        # Plot residulas vs fitted values
        # fig_res, ax_res = plt.subplots()
        ax_res = ax_all[0,1]
        ax_res.scatter(fitted_values, residuals)
        ax_res.hlines(y=0, xmin=min(fitted_values), xmax=max(fitted_values), linestyle = 'dashed',\
                      color = 'black')
        ax_res.set_xlabel('Fitted Values')
        ax_res.set_ylabel('Residuals')
        ax_res.set_title('Versus Fit')
        # ax_res.set_xlim(min(fitted_values)-0.0001,max(fitted_values)+0.0001)
        
        # Plot histogram of residuals
        # fig_hist, ax_hist = plt.subplots()
        ax_hist = ax_all[1,0]
        ax_hist.hist(residuals, color='#0504aa', alpha=0.7, rwidth=0.9)
        ax_hist.set_xlabel('Residuals')
        ax_hist.set_ylabel('Counts')
        ax_hist.set_title('Histogram')

        # Plot residuals vs order
        # fig_order, ax_order = plt.subplots()
        ax_order = ax_all[1,1]
        order_list = np.arange(1,len(residuals)+1,1)
        ax_order.plot(order_list,residuals,'-o')
        ax_order.hlines(y=0, xmin=min(order_list), xmax=max(order_list), linestyle = 'dashed',\
                      color = 'black')
        ax_order.set_xlabel('Observation Order')
        ax_order.set_ylabel('Residuals')
        ax_order.set_title('Versus Order')

    plt.show()


if __name__ == '__main__':
    alpha_values = list(df_fromexcel['alpha'])
    k_max = 30
    min_t_ent_plot(alpha_values, k_max)
    Me_t_ent_plot(alpha_values, k_max)
    mean_t_ent_plot(alpha_values, k_max)

