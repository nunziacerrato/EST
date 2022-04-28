''' In this program the theoretical probability distribution functions that have been found with the
    Minitab software (the 3-parameter Gamma distribution and the 3-parameter Lognormal distribution)
    are overlapped in the same graph, together with the data. The parameter reported for the three
    cases considered here are the ones obtained in the analysis with Minitab.
'''
import numpy as np
import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt

# Parameter
N = 2
iterations = 20000
fixed, fixed_val = 'g', 1
gamma = fixed_val
common_path = "C:\\Users\\cerra\\Documents\\GitHub\\lindblad_data"

# Individual plot for different k values with both theoretical pdfs overlapped
if False:
    k_list = [0,1,1000]

    # Lognormal parameter
    scale_location_lognormal = [-1.32915,-1.86269,-3.35909] # scale = exp(location)
    loc_threshold_lognormal = [0.78692,0.80338,0.82252] # loc
    sigma_scale_lognormal = [0.49121,0.54827,1.14009] # sigma or s

    # Gamma parameter
    shape_gamma = [2.89677,2.45291,0.89219]
    loc_threshold_gamma = [0.83435,0.82813,0.82400]
    scale_gamma = [0.08656,0.06334,0.06776]

    for k_value in k_list:
        alpha = k_value*gamma
        fig_hist, ax_hist = plt.subplots(figsize=(15,10))
        ax_hist.set_title(fr'$\tau_{{ent}}$ distribution for k={k_value} with theoretical fit pdfs',\
                            fontsize = 18)
        # Create a DataFrame from the excel file where all the data are saved
        df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')

        # Access to the column 't_ent' and drop NaN values if present
        t_ent_values = df_fromexcel['t_ent'].dropna().values
        t_ent_range = np.linspace(min(t_ent_values), max(t_ent_values), 500)

        # Compute the histogram of t_ent values
        ax_hist.hist(t_ent_values, bins = 'auto', histtype='bar', ec="black", fc="cornflowerblue",\
            alpha=0.5, density = True, label = fr'Data - k = {k_value}, $\alpha$ = {np.round(alpha,2)},'\
            fr' $\gamma$ = {np.round(gamma,2)}')
        
        # Access to the Lognormal parameter and plot the Lognormal pdf
        sigma = sigma_scale_lognormal[k_list.index(k_value)]
        loc = loc_threshold_lognormal[k_list.index(k_value)]
        scale = np.exp(scale_location_lognormal[k_list.index(k_value)])

        rv_lognorm = scipy.stats.lognorm(sigma, loc = loc, scale = scale)
        ax_hist.plot(t_ent_range, rv_lognorm.pdf(t_ent_range), 'blue', linewidth=3,\
                    label='3-parameter Lognormal pdf')

        # Access to the Gamma parameter and plot the Gamma pdf
        shape = shape_gamma[k_list.index(k_value)]
        loc = loc_threshold_gamma[k_list.index(k_value)]
        scale = scale_gamma[k_list.index(k_value)]

        rv_gamma = scipy.stats.gamma(shape, loc = loc, scale = scale)
        if k_value == 1000:
            t_ent_range = np.linspace(min(t_ent_values)+0.0015, max(t_ent_values), 500)
        # scipy.stats.gamma.ppf(0.01, shape, loc = loc, scale = scale)
        ax_hist.plot(t_ent_range, rv_gamma.pdf(t_ent_range), 'red', linewidth=3,\
                    label='3-parameter Gamma pdf')    

        
        ax_hist.set_xlabel(fr'$\tau_{{ent}}$', fontsize = 16)
        ax_hist.set_ylabel('Density', fontsize = 14)
        ax_hist.legend(fontsize = 14)
        plt.show()

# Layout 3x2 with both theoretical pdfs overlapped
if True:
    k_list = [0,0.5,1,10,100,1000]
    fig_hist, ax_hist = plt.subplots(nrows=3, ncols=2)
    
    # Lognormal parameter
    scale_location_lognormal = [-1.32915,-1.5426,-1.86269,-3.30639,-3.36680,-3.35909] # scale = exp(location)
    loc_threshold_lognormal = [0.78692,0.79677,0.80338,0.82249,0.82264,0.82252] # loc
    sigma_scale_lognormal = [0.49121,0.49246,0.54827,1.10641,1.15127,1.14009] # sigma or s

    # Gamma parameter
    shape_gamma = [2.89677,2.89435,2.45291,0.93225,0.88896,0.89219]
    loc_threshold_gamma = [0.83435,0.83462,0.82410,0.82813,0.82400,0.82400]
    scale_gamma = [0.08656,0.07022,0.06334,0.06660,0.06836,0.06776]

    for k_value in k_list:
        alpha = k_value*gamma
        k_index = k_list.index(k_value)
        
        fig_hist.suptitle(fr'Layout $\tau_{{ent}}$ distributions with theoretical fit pdfs overlapped',\
                            fontsize = 18)
        # Create a DataFrame from the excel file where all the data are saved
        df_fromexcel = pd.read_excel(f'{common_path}\\{fixed}_fixed_comparable_data\\{fixed}={fixed_val}\\'\
            f'Lindblad_eigvals_t_ent_{fixed}_fixed_k_{k_value}_{iterations}_iterations_N_{N}.xlsx')

        # Access to the column 't_ent' and drop NaN values if present
        t_ent_values = df_fromexcel['t_ent'].dropna().values
        t_ent_range = np.linspace(min(t_ent_values), max(t_ent_values), 500)

        # Compute the histogram of t_ent values
        ax_hist[k_index//2, k_index%2].hist(t_ent_values, bins = 'auto', histtype='bar', ec="black",\
            fc="cornflowerblue", alpha=0.5, density = True,\
            label = fr'Data')# - k = {k_value}, $\alpha$ = {np.round(alpha,2)},'\
            #fr' $\gamma$ = {np.round(gamma,2)}')
        
        # Access to the Lognormal parameter and plot the Lognormal pdf
        sigma = sigma_scale_lognormal[k_index]
        loc = loc_threshold_lognormal[k_index]
        scale = np.exp(scale_location_lognormal[k_index])

        rv_lognorm = scipy.stats.lognorm(sigma, loc = loc, scale = scale)
        ax_hist[k_index//2, k_index%2].plot(t_ent_range, rv_lognorm.pdf(t_ent_range), 'blue',\
                                            linewidth=2, label='3-parameter Lognormal pdf')

        # Access to the Gamma parameter and plot the Gamma pdf
        shape = shape_gamma[k_index]
        loc = loc_threshold_gamma[k_index]
        scale = scale_gamma[k_index]

        rv_gamma = scipy.stats.gamma(shape, loc = loc, scale = scale)
        if k_value == 1000 or k_value == 100:
            t_ent_range = np.linspace(min(t_ent_values)+0.0015, max(t_ent_values), 500)
        # scipy.stats.gamma.ppf(0.01, shape, loc = loc, scale = scale)
        ax_hist[k_index//2, k_index%2].plot(t_ent_range, rv_gamma.pdf(t_ent_range), 'red',\
                                            linewidth=2, label='3-parameter Gamma pdf')    

        ax_hist[k_index//2, k_index%2].set_title(fr'k={k_value}, $\alpha$={alpha}, $\gamma$={gamma}',\
                                                fontsize = 14)
        ax_hist[k_index//2, k_index%2].set_xlabel(fr'$\tau_{{ent}}$', fontsize = 13)
        ax_hist[k_index//2, k_index%2].set_ylabel('Density', fontsize = 12)
        ax_hist[k_index//2, k_index%2].legend(fontsize = 10)
        for ax in ax_hist.flat:
            ax.label_outer()
        # ax_hist.legend(fontsize = 14)
    plt.show()

