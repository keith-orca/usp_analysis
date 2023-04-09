import scipy.stats as stats
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from usp_functions import *
from scipy.stats import pearsonr
from usp_file_writing import to_excel_named


def combinations_anova(vpp, f_mod, photometry_data):
    # sorts 2D photometry arrays by indexed vpp f_mod combination and averages within the sorted groups.
    # This can be used on 2D (photometry fluorescence data) or 1D arrays (slice mean arrays)
    try:
        combinations, vpp_set, f_mod_set = get_parameter_combinations(vpp, f_mod)
        combination_photometry_data = [[] for i in range(len(combinations))]
        average_set = []

        for i in range(len(photometry_data)):
            combination_index = combinations.index([vpp[i], f_mod[i]])
            combination_photometry_data[combination_index].append(photometry_data[i])

        statistic, p_val = stats.f_oneway(*combination_photometry_data)
        print('p= ' + str(p_val))
        print(statistic)

        plt.figure()
        sns.swarmplot(data=combination_photometry_data[:], color="0", alpha=.35)
        sns.violinplot(data=combination_photometry_data[:], capsize=.1, ci="sd")

        for i in range(len(combination_photometry_data)):
            average_set.append(sum(combination_photometry_data[i]) / len(combination_photometry_data[i]))

        avg_set = [f_mod_set, average_set]
        to_excel_named(avg_set, 'avg_set')
    except:
        print('Insufficient group n for statistical analysis.')


def pearson_correlate(list_1,list_2):
    corr, p_value = pearsonr(list_1, list_2)
    r_squared = corr*corr

    return r_squared, p_value