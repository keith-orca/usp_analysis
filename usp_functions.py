# usp analysis functions
import numpy as np
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
import copy
import usp_file_writing
from scipy.optimize import curve_fit


def get_parameter_combinations(vpp, f_mod):
    # Iterates through vpp and f_mod, appending any combination of the two not already in present in
    # the combinations list. Unique combinations are appended as a list [vpp, f_mod]

    combinations = []

    for i in range(len(vpp)):
        combination = [vpp[i], f_mod[i]]
        if not combination in combinations:
            combinations.append(combination)

    combinations = sorted(combinations)
    vpp_set = [combinations[i][0] for i in range(len(combinations))]
    f_mod_set = [combinations[i][1] for i in range(len(combinations))]

    return combinations, vpp_set, f_mod_set


def get_combination_means(vpp, f_mod, photometry_data, output_name=False, write_raw_combinations=False):
    # sorts 2D photometry arrays by indexed vpp f_mod combination and averages within the sorted groups.
    # This can be used on 2D (photometry fluorescence data) or 1D arrays (slice mean arrays)

    combinations, vpp_set, f_mod_set = get_parameter_combinations(vpp, f_mod)
    combination_photometry_data = [[] for i in range(len(combinations))]

    for i in range(len(photometry_data)):
        combination_index = combinations.index([vpp[i], f_mod[i]])
        combination_photometry_data[combination_index].append(photometry_data[i])

    combination_means = []
    combination_ses = []

    for i in range(len(combination_photometry_data)):
        combination_means.append(np.mean(combination_photometry_data[i], axis=0))
        combination_ses.append(stats.sem(combination_photometry_data[i], axis=0))

    if output_name:

        for i in range(len(combinations)):
            combinations[i].append(combination_means[i])

        usp_file_writing.to_excel_named(combinations, output_name)

    if write_raw_combinations:
        for comb_index in range(len(combination_photometry_data)):
            comb_name = str(combinations[comb_index][0]) + "_" + str(combinations[comb_index][1])
            to_excel_named(transpose(combination_photometry_data[comb_index]), comb_name)

    return combination_means, combination_ses


def get_collective_mean(photometry_data):
    # sorts 2D photometry arrays by indexed vpp f_mod combination and averages within the sorted groups.
    # This can be used on 2D (photometry fluorescence data) or 1D arrays (slice mean arrays)

    mean = np.mean(photometry_data, axis=0)
    sem = stats.sem(photometry_data, axis=0)
    mean = [mean]
    sem = [sem]

    return mean, sem


def bin_vector(vector, bin_size, numpy_array=False):

    binned_vector = []
    for i in range(0, len(vector),bin_size):
        binned_vector.append(sum(vector[i:int(i+bin_size)])/bin_size)

    if numpy_array:
        binned_vector = np.asarray(binned_vector, dtype=np.float32)

    return binned_vector


def calculate_delta_f(photometry_data, baseline_start, baseline_end, fs):
    # Calculates delta f for all trials of photometry data using the given baseline slice for normalization

    delta_f = np.copy(photometry_data)
    empty_trial_count = 0

    for i in range(len(delta_f)):
        baseline_mean = delta_f[i][int(baseline_start * fs): int(baseline_end * fs)]
        if sum(baseline_mean) != 0:
            baseline_mean = baseline_mean.mean()
            delta_f[i] = ((delta_f[i] - baseline_mean) / baseline_mean) * 100
        else:
            empty_trial_count = empty_trial_count + 1
            delta_f[i] = 0
    if empty_trial_count != 0:
        print('Error: Found ' + str(empty_trial_count) + ' trial/s with empty data set.')

    return delta_f


def calculate_delta_f(photometry_data, baseline_start, baseline_end, fs):
    # Calculates delta f for all trials of photometry data using the given baseline slice for normalization

    delta_f = np.copy(photometry_data)
    empty_trial_count = 0

    for i in range(len(delta_f)):
        baseline_mean = delta_f[i][int(baseline_start * fs): int(baseline_end * fs)]
        if sum(baseline_mean) != 0:
            baseline_mean = baseline_mean.mean()
            delta_f[i] = ((delta_f[i] - baseline_mean) / baseline_mean) * 100
        else:
            empty_trial_count = empty_trial_count + 1
            delta_f[i] = 0
    if empty_trial_count != 0:
        print('Error: Found ' + str(empty_trial_count) + ' trial/s with empty data set.')

    return delta_f


def strip_file_path(file_path_name):
    # function for removing folder names from file_path_name by taking all chars from end to last '/'

    stripped_file_name = ''

    for i in range(1, len(file_path_name) - 1):
        if file_path_name[-i] != '/':
            stripped_file_name = file_path_name[-i] + stripped_file_name
        else:
            break

    return stripped_file_name


def get_whole_trace(usp_instance):
    trace_set = []
    for i in range(len(usp_instance.photometry_ch)):
        trace_set.append([])
        for n in range(len(usp_instance.photometry_ch[i])):
            trace_set[i].extend(usp_instance.photometry_ch[i][n])

    return trace_set


def moving_average(values, window):
    values_padded = np.pad(values, (window // 2, window - 1 - window // 2), mode='edge')
    weights = np.repeat(1.0, window) / window
    values_smooth = np.convolve(values_padded, weights, mode='valid')

    return values_smooth


def exp_decay_func(x, a, k, b):
    return a * np.exp(-k*x) + b


def lin_decay_func(x, m, b):
    return m * x + b


def remove_exp_decay(data):
    # assumes linear x and 2D data of trials from continuous recording

    x = np.linspace(0, len(data), len(data))
    y = data
    p0 = (1, 0, 1) # starting search koefs
    opt, pcov = curve_fit(exp_decay_func, x, y, p0, maxfev=1200)
    a, k, b = opt
    # test result
    fit_x = np.linspace(0, len(data), len(data))
    fit_y = exp_decay_func(fit_x, a, k, b)

    return fit_y


def remove_lin_decay(data):
    # assumes linear x and 2D data of trials from continuous recording

    x = np.linspace(0, len(data), len(data))
    y = data

    m, b = np.polyfit(x, y, 1)

    fit_x = np.linspace(0, len(data), len(data))
    fit_y = lin_decay_func(fit_x, m, b)

    return fit_y


def transpose(matrix):
    # Flips matrix rows into columns and vice versa. Matrix should contain lists of = length,
    # otherwise short lists will be extended with ""
    listLength = []
    for i in range(len(matrix)):
        listLength.append(len(matrix[i]))
    long = max(listLength)
    for i in range(len(matrix)):
        while len(matrix[i]) < long:
            matrix[i].append("")
    longList = []
    for i in range(len(matrix[0])):
        shortList = []
        for n in range(len(matrix)):
            shortList.append(matrix[n][i])
        longList.append(shortList)
    return longList


def get_abs_max_2d(list_2d):
    max_abs_val = False
    for i in range(len(list_2d)):
        for n in range(len(list_2d[i])):
            if abs(list_2d[i][n]) > max_abs_val:
                max_abs_val = abs(list_2d[i][n])

    return max_abs_val


def get_abs_mean_2d(list_2d):
    means_list = []
    for i in range(len(list_2d)):
        for n in range(len(list_2d[i])):
            means_list.append(np.mean(abs(list_2d[i][n])))
    mean_abs_val = np.mean(means_list)

    return mean_abs_val


def print_parameters(usp_instance):
    print("t_pre: " + str(usp_instance.t_pre))
    print("t_on: " + str(usp_instance.t_on))
    print("t_post: " + str(usp_instance.t_post))
    print("f_c: " + str(usp_instance.f_c))
    print("fs: " + str(usp_instance.fs))
    print("duty_cycle: " + str(usp_instance.duty_cycle))
    print("animal id: " + str(usp_instance.animal_id))

