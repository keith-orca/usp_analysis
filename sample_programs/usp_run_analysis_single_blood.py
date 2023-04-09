import usp_load_data as lmd
from usp_graphics import *
from usp_signal_analysis import *
import usp_statistics
from usp_file_writing import *

# load several .mat files from the us_photometry matlab data acquisition software
#usp = lmd.load_folder('C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/smst_641')

# channel of interest for plotting; channel 0 is blue channel or photometry_ch1 from matlab file
ch = 0

# can also load .npy files written by the usp_photometry class method write_file  (doesn't get used in this example)

usp = lmd.load_npy_file('C:\\Users\\15163\\Desktop\\stanford_software\\usp_analysis_code\\usp_data\\bns6_prf_1.npy')
#usp2 = lmd.load_npy_file('C:\\Users\\15163\\Desktop\\stanford_software\\usp_analysis_code\\usp_data\\LC765_durationintensity_1.npy')
#usp.acquire_usp_data(usp2)
#usp.write_file('cmtb14_20hz40s_intensity_new')

# turn ON for blood expts
usp.remove_decay(ch, padding=True, decay_type='exponential') # padding offsets - values by adding back the abs(min) of the whole trace

usp.bin_data(4)
#usp.remove_events(ch, 1)
#usp.subtract_out(0.2)

#usp.filter_data(filter_vals=[0.5], must_contain_all=False)
c, vpp_set, f_mod_set = get_parameter_combinations(usp.vpp, usp.f_mod)
trace = get_whole_trace(usp)
line_plot(trace[ch])
print_parameters(usp)

# smooth channel data with moving average
#usp.smooth_signal(0, window=5)

# converts original and filtered data (if existing) to delta baseline where baseline is within input range in seconds
usp.get_delta_f(0, usp.t_pre)

# raise the gain on the UV channel since it does not fully correct for gap loss at present levels
#usp.adjust_gain_delta_f(1, gain=1)

# Get a single list of arrays (note [ch] does not exist for usp.subtracted_delta_f
# for blood, no normalizing and no delta subtration, but added in exp decay subtraction, switch subtract_delta order if used
#usp.subtract_delta_f(0, 1)
usp.normalize_delta_f(channel=ch, norm_base='mean')

#usp.get_fft(ch, 0, 5)
#fft_combination_means, fft_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.fft_vals)
#shade_means_plot(fft_combination_means, fft_combination_ses, usp, usp.vpp, usp.f_mod, x_axis=usp.fft_freqs)
#usp.remove_events_by_phase(0, 1, criteria=0)

# get the means of trials within photometry_ch1 using the combinations of vpp and f_mod present
combination_means, combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.photometry_ch[ch])

# do the same thing for the delta transforms on photometry_ch
delta_f_combination_means, delta_f_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f[ch])

# Get means of all photometry trials within indicated window, take means of means categorized by vpp-f_mod combination
# and then plot the vpp values with corresponding photometry means.

#usp.get_slice_means(usp.t_pre, usp.t_pre+usp.t_on)
usp.get_slice_means(5, 10)
slice_combination_means, slice_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f_slice_means[ch], output_name='output')
combinations, vpp_set, f_mod_set = get_parameter_combinations(usp.vpp, usp.f_mod)
plot_xy(vpp_set, slice_combination_means, x_axis_title='voltage')
to_excel_named(transpose(usp.delta_f_slice_means), 'slicemeans')
usp_statistics.combinations_anova(usp.vpp, usp.f_mod, usp.delta_f_slice_means[ch])

#interpolate_grid(combinations, slice_combination_means)

#avg_lines_stacked(delta_f_combination_means, usp, usp.vpp, usp.f_mod)

#means, sess = get_collective_mean(usp.vpp_filtered, usp.f_mod_filtered, usp.delta_f_filtered[ch])

shade_means_plot(delta_f_combination_means, delta_f_combination_ses, usp, usp.vpp, usp.f_mod)
to_excel_named(transpose(delta_f_combination_means), 'deltaf_combination_means')
all_lines_stacked(delta_f_combination_means, usp, usp.vpp, usp.f_mod, channel=ch, delta=True)

# Get a visual of stimulation series across trials
#stimulation_map(usp)

# print the file names composing the input usp object
print(usp.file_names)

# Show all the figures generated within the functions called from usp_graphics.py
plt.show()

# Write the usp object to a numpy dictionary file; this is useful if files have been merged
#usp.write_file('cmtb14_20hz40s_intensity_combine')
