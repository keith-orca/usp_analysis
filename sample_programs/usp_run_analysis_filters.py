import usp_load_data as lmd
from usp_graphics import *
import usp_statistics

# load several .mat files from the us_photometry matlab data acquisition software

# can also load .npy files written by the usp_photometry class method write_file  (doesn't get used in this example)
#usp = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_194carrier_2014_8_22_13_29.mat')
#usp2 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_494carrier_2014_8_22_14_13.mat')
#usp3 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_804carrier_2014_8_22_15_5.mat')
#usp4 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_1600000_carrier_2014_8_22_15_59.mat')
#usp5 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_1230000_carrier_2014_8_23_9_27.mat')
#usp6 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_1820000_carrier_2014_8_22_18_1.mat')
#usp7 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_804000_carrier_run2_2014_8_23_10_45.mat')
#usp8 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cst#13_gcamp6s_cont_wave_2s_1072000_carrier_2014_8_23_11_30.mat')

usp = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/ca1_camkII_gcamp7s_np_jordan_cont_wave_2s_2014_8_22_8_29.mat')
usp2 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/ca1_camkII_gcamp7s_cont_wave_2s_1234000_carrier_2014_8_23_12_40.mat')
usp3 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/ca1_camkII_gcamp7s_cont_wave_2s_1820000_carrier_2014_8_23_15_41.mat')
usp4 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/ca1_camkII_gcamp7s_cont_wave_2s_1600000_run2_carrier_2014_8_23_16_54.mat')
usp5 = lmd.load_mat_data('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/ca1_camkII_gcamp7s_cont_wave_2s_194000_carrier_2014_8_23_18_6.mat')

# This method merges data from one usp object into another. Trials from each file are still separable by
# looking at the file_index attribute of the object where the index points to the associated file_name
usp.acquire_usp_data(usp2)
usp.acquire_usp_data(usp3)
usp.acquire_usp_data(usp4)
usp.acquire_usp_data(usp5)


# Get a copy of data in which events are filtered to include events that have at least one parameter from the input list
usp.filter_data(filter_vals=[1820000,0.0001,0.2, 0.4, 0.8, 0.6,0.8], must_contain_all=True)

# converts original and filtered data (if existing) to delta baseline where baseline is within input range in seconds
usp.get_delta_f(0, 1)

c, vpp_set, f_mod_set = get_parameter_combinations(usp.vpp, usp.f_mod)
print(vpp_set)
print(f_mod_set)

# get the means of trials within photometry_ch1 using the combinations of vpp and f_mod present
combination_means, combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.photometry_ch1)

# do the same thing for the delta transforms on photometry_ch1
delta_f_combination_means, delta_f_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f)

# do the same thing for the filtered delta transforms on photometry_ch1
delta_f_filtered_combination_means, delta_f_filtered_combination_ses = get_combination_means(usp.vpp_filtered,
                                                                                             usp.f_mod_filtered,
                                                                                             usp.delta_f_filtered)

# Get means of all photometry trials within indicated window, take means of means categorized by vpp-f_mod combination
# and then plot the vpp values with corresponding photometry means.
usp.get_slice_means(1, 3)
slice_combination_means, slice_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f_slice_means)
combinations, vpp_set, f_mod_set = get_parameter_combinations(usp.vpp, usp.f_mod)
plot_xy(vpp_set, slice_combination_means, x_axis_title='voltage')

usp_statistics.combinations_anova(usp.vpp, usp.f_mod, usp.delta_f_slice_means)

interpolate_grid(combinations, slice_combination_means)

#avg_lines_stacked(delta_f_filtered_combination_means, usp, usp.vpp, usp.f_mod)

shade_means_plot(delta_f_filtered_combination_means, delta_f_filtered_combination_ses, usp, usp.vpp_filtered,
                 usp.f_mod_filtered)

shade_means_plot(delta_f_combination_means, delta_f_combination_ses, usp, usp.vpp,
                 usp.f_mod)

all_lines_stacked(delta_f_filtered_combination_means, usp, usp.vpp_filtered, usp.f_mod_filtered, delta=True,
                  filtered=True)

# Get a visual of stimulation series across trials
stimulation_map(usp)

# print the file names composing the input usp object
print(usp.file_names)

# Show all the figures generated within the functions called from usp_graphics.py
plt.show()

# Write the usp object to a numpy dictionary file; this is useful if files have been merged
usp.write_file('hcrt_expt_1_2')
