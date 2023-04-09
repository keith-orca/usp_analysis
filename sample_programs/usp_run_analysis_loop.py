import usp_load_data as lmd
from usp_graphics import *

file_id = ['C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cort-gcamp6f-2-longrun_2014_8_2_19_7.mat',
           'C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cort-gcamp6f-2-longrun_2014_8_2_18_41.mat',
           'C:/Users/Keith/Desktop/usp_analysis_code/usp_sample_data/cort-gcamp6f-2-longrun_2014_8_2_17_55.mat']
           
           
           
for i in range(len(file_id)):
    # load several .mat files from the us_photometry matlab data acquisition software
    usp = lmd.load_mat_data(file_id[i])


    # Get a copy of data in which events are filtered to include events that have at least one parameter from the input list
    usp.filter_data(filter_vals=[240000])

    # converts original and filtered data (if existing) to delta baseline where baseline is within input range in seconds
    usp.get_delta_f(0, 1)

    # get the means of trials within photometry_ch1 using the combinations of vpp and f_mod present
    combination_means, combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.photometry_ch1)

    # do the same thing for the delta transforms on photometry_ch1
    delta_f_combination_means, delta_f_combination_ses = get_combination_means(usp.vpp, usp.f_mod, usp.delta_f)

    # do the same thing for the filtered delta transforms on photometry_ch1
    delta_f_filtered_combination_means, delta_f_filtered_combination_ses = get_combination_means(usp.vpp_filtered,
                                                                                                 usp.f_mod_filtered,
                                                                                                 usp.delta_f_filtered)


    avg_lines_stacked(delta_f_filtered_combination_means, usp, usp.vpp, usp.f_mod)

    shade_means_plot(delta_f_filtered_combination_means, delta_f_filtered_combination_ses, usp, usp.vpp_filtered,
                     usp.f_mod_filtered)

    all_lines_stacked(delta_f_filtered_combination_means, usp, usp.vpp_filtered, usp.f_mod_filtered, delta=True,
                      filtered=True)

    # Get a visual of stimulation series across trials
    stimulation_map(usp)

    # print the file names composing the input usp object
    print(usp.file_names)

    # Show all the figures generated within the functions called from usp_graphics.py
    plt.show()

