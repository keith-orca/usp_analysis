import usp_load_data as lmd
from usp_graphics import *

# load several .mat files from the us_photometry matlab data acquisition software
usp = lmd.load_npy_file('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cmtb17_2.5Hz_5s_intensity_1.npy')


# can also load .npy files written by the usp_photometry class method write_file  (doesn't get used in this example)
#usp3 = lmd.load_npy_file('C:/Users/Keith/Desktop/usp_analysis_code/f_test.npy')

# This method merges data from one usp object into another. Trials from each file are still separable by
# looking at the file_index attribute of the object where the index points to the associated file_name
#usp.acquire_usp_data(usp2)

# Get a copy of data in which events are filtered to include events that have at least one parameter from the input list
#usp.filter_data(filter_vals=[10000])

# converts original and filtered data (if existing) to delta baseline where baseline is within input range in seconds
usp.get_delta_f(0, 1)

# get the means of trials within photometry_ch1 using the combinations of vpp and f_mod present
combination_means, combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.photometry_ch1)

# do the same thing for the delta transforms on photometry_ch1
delta_f_combination_means, delta_f_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f)

# do the same thing for the filtered delta transforms on photometry_ch1
#delta_f_filtered_combination_means, delta_f_filtered_combination_ses = get_combination_means(usp.vpp_filtered,
#                                                                                             usp.f_mod_filtered,
#                                                                                             usp.delta_f_filtered)

# Get means of all photometry trials within indicated window, take means of means categorized by vpp-f_mod combination
# and then plot the vpp values with corresponding photometry means.
slice_means = get_slice_means(usp.delta_f, 2, 2.3, usp.fs)
slice_combination_means, slice_combination_ses = get_combination_means(usp.vpp, usp.f_mod,
                                                                       slice_means)
combinations, vpp_set, f_mod_set = get_parameter_combinations(usp.vpp, usp.f_mod)
#plot_xy(vpp_set, slice_combination_means, x_axis_title='voltage')

interpolate_grid(combinations, slice_combination_means)

#avg_lines_stacked(delta_f_filtered_combination_means, usp, usp.vpp, usp.f_mod)

shade_means_plot(delta_f_combination_means, delta_f_combination_ses, usp, usp.vpp,
                 usp.f_mod)

all_lines_stacked(delta_f_combination_means, usp, usp.vpp, usp.f_mod, delta=True,
                  filtered=False)

# Get a visual of stimulation series across trials
stimulation_map(usp)

# print the file names composing the input usp object
print(usp.file_names)

# Show all the figures generated within the functions called from usp_graphics.py
plt.show()

# Write the usp object to a numpy dictionary file; this is useful if files have been merged
usp.write_file('f_test')
