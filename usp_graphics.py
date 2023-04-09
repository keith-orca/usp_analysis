# Module for plotting functions
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib import cm
from numpy import mgrid
from usp_functions import *


def shade_means_plot(combination_means, combination_ses, us_photometry_instance, vpp, f_mod, x_axis=False, y_axis_scale=False):
    # Create shaded error bars plot
    plt.figure()
    combinations, vpp_set, f_mod_set = get_parameter_combinations(vpp, f_mod)

    if x_axis == False:
        time = np.arange(0, len(combination_means[0]) / us_photometry_instance.fs, 1 / us_photometry_instance.fs)
    else:
        time = x_axis[0]

    print(combinations)
    for i in range(len(combination_means)):
        plt.plot(time, combination_means[i])
    plt.legend(combinations, framealpha=0, loc='upper left', fontsize=7)

    for i in range(len(combination_means)):
        plt.fill_between(time, combination_means[i] - combination_ses[i], combination_means[i] + combination_ses[i],
                         alpha=0.4)
    ax = plt.gca()
    ax.set_ylabel('ΔF/F', fontsize=13, color=[0.25, 0.25, 0.25])
    ax.set_xlabel('Time (s)', fontsize=13, color=[0.25, 0.25, 0.25])
    ax.tick_params(axis='both', which='major', width=1.1, labelsize=10, labelcolor=[0.3, 0.3, 0.3],
                   color=[0.4, 0.4, 0.4])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for spine in ax.spines.values():
        spine.set_edgecolor([0.4, 0.4, 0.4])
        spine.set_linewidth(1.1)

    if y_axis_scale != False:
        plt.yscale(y_axis_scale)

    ax = plt.gca()
    #ax.set_ylabel('ΔF/F')
    #ax.set_xlim(0, 10)
    #ax.set_ylim(-1, 1.5)


def avg_lines_stacked(combination_means, us_photometry_instance, vpp, f_mod):
    # Create stack of average line plots

    plt.figure()
    combinations, vpp_set, f_mod_set = get_parameter_combinations(vpp, f_mod)
    time = np.arange(0, len(combination_means[0]) / us_photometry_instance.fs, 1 / us_photometry_instance.fs)

    for i in range(len(combination_means)):
        plt_nums = int(str(len(combination_means)) + '1' + str(i + 1))
        plt.subplot(plt_nums)
        plt.plot(time, combination_means[i])
        plt.title(str(combinations[i]))
        ax = plt.gca()
        ax.set_ylabel('ΔF/F')

    ax.set_xlabel('Time (s)')
    plt.subplots_adjust(hspace=0.7)


def all_lines_stacked(combination_means, us_photometry_instance, vpp, f_mod, channel=0, delta=False):
    # Create plot with all trials and mean (colored line) for each combination where each combination has its
    # own subplot

    plt.figure()
    combinations, vpp_set, f_mod_set = get_parameter_combinations(vpp, f_mod)
    time = np.arange(0, len(combination_means[0]) / us_photometry_instance.fs, 1 / us_photometry_instance.fs)

    for i in range(len(combination_means)):
        plt.subplot(int(len(combination_means)), 1, (i + 1))

        if delta:
            lines = np.copy(us_photometry_instance.delta_f[channel])
        else:
            lines = np.copy(us_photometry_instance.photometry_ch[channel])
        color_map = cm.get_cmap('inferno', len(lines))

        for n in range(len(lines)):
            if vpp[n] == combinations[i][0] and f_mod[n] == combinations[i][1]:
                #plt.plot(time, lines[n], color=color_map(n), alpha=1, linewidth=0.3)
                plt.plot(time, lines[n], color='black', alpha=0.3, linewidth=0.7)

        plt.plot(time, combination_means[i], color=[.7, 0, 0.5], alpha=0.5, linewidth=2)
        plt.title(str(combinations[i]), fontsize=10)
        ax = plt.gca()
        ax.set_ylabel('ΔF/F')
        #ax.set_xlim(0, 20)
        #ax.set_ylim(-5, 11)

    ax.set_xlabel('Time (s)')  # gives x axis to last plot in loop at bottom
    plt.subplots_adjust(hspace=0.7)


def plot_xy(x, y, x_axis_title='vpp or f_mod', y_axis_title='Fluoresence (mV)'):
    plt.figure()
    ax = plt.gca()
    ax.set_ylabel(y_axis_title)
    ax.set_xlabel(x_axis_title)
    ax = plt.gca()
    #ax.set_ylabel('ΔF/F')
    #ax.set_xlim(0, 10)
    #ax.set_ylim(-0.5, 2)
    plt.plot(x, y)


def stimulation_map(us_photometry_instance):
    # Creates a 2d array of vpp and f_mod for time points within each trial for all trials
    # and displays data with imshow()

    usp = us_photometry_instance
    vpp_map = np.zeros((len(usp.vpp), int((usp.t_pre + usp.t_on + usp.t_post) * usp.fs)))
    f_mod_map = np.zeros((len(usp.f_mod), int((usp.t_pre + usp.t_on + usp.t_post) * usp.fs)))
    time = np.arange(0, len(usp.vpp) / usp.fs, 1 / usp.fs)

    for i in range(len(usp.vpp)):
        stim_pre = np.zeros(int(usp.t_pre * usp.fs))
        stim_on_vpp = np.full((int(usp.t_on * usp.fs)), usp.vpp[i])
        stim_on_f_mod = np.full((int(usp.t_on * usp.fs)), usp.f_mod[i])
        stim_post = np.zeros(int(usp.t_post * usp.fs))
        vpp_part_map = np.append(stim_pre, stim_on_vpp)
        f_mod_part_map = np.append(stim_pre, stim_on_f_mod)
        vpp_map[i] = np.append(vpp_part_map, stim_post)
        f_mod_map[i] = np.append(f_mod_part_map, stim_post)

    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(9, 4), sharey=True)
    ax1.set_title('Voltage')
    ax2.set_title('Frequency modulation')
    ax1.set_ylabel('Trial #')
    ax1.set_xlabel('Time')
    ax2.set_xlabel('Time')
    vpp_ax = ax1.imshow(vpp_map, cmap='gist_gray_r')
    f_mod_ax = ax2.imshow(f_mod_map, cmap='gist_gray_r')
    fig.colorbar(vpp_ax, ax=ax1)
    fig.colorbar(f_mod_ax, ax=ax2)


def make_cmap():
    # just sticking this here temporarily
    top = cm.get_cmap('Oranges_r', 128)
    bottom = cm.get_cmap('Blues', 128)

    newcolors = np.vstack((top(np.linspace(0, 1, 128)),
                           bottom(np.linspace(0, 1, 128))))
    newcmp = ListedColormap(newcolors, name='OrangeBlue')


def interpolate_grid(points, values):
    # Takes points (array of [x,y] arrays), and values for each point (array), and interpolates onto
    # all positions onto a point map given by np.mgrid. May need at least ~ 6 data points

    try:
        points = np.array(points)
        vpp_min = points.transpose()[0].min()
        vpp_max = points.transpose()[0].max()
        f_mod_min = points.transpose()[1].min()
        f_mod_max = points.transpose()[1].max()
        grid_x, grid_y = np.meshgrid(np.linspace(vpp_min, vpp_max, 1000), np.linspace(f_mod_min, f_mod_max, 1000))
        grid_z = griddata(points, values, (grid_x, grid_y), method='cubic')
        plt.figure()
        plt.imshow(grid_z, extent=(vpp_min, vpp_max, f_mod_min, f_mod_max), aspect='auto', origin='lower',
                   cmap='inferno')
        ax = plt.gca()
        ax.ticklabel_format(axis='both', style='sci', scilimits=(-2, 2))
        ax.set_ylabel('f_mod')
        ax.set_xlabel('vpp')
        plt.colorbar()
        plt.scatter(points.transpose()[0], points.transpose()[1], c='black', s=20, alpha=0.3, marker="o",
                    edgecolors='none')

    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)


def histogram(data_set):
    plt.figure()
    plt.hist(data_set, bins=100)


def line_plot(vector):
    plt.figure()
    plt.plot(vector, color='black', alpha=1, linewidth=0.1)


def scatter_plot(x,y):
    plt.figure()
    plt.scatter(x, y, alpha=0.5)