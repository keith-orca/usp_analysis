from usp_functions import *
from usp_signal_analysis import *
import numpy as np
import copy
import time
import matplotlib.pyplot as plt
import usp_statistics
from usp_graphics import *
from usp_file_writing import *

class USPhotometry:

    def __init__(self, t_pre, t_on, t_post, fs, f_c, f_mod, vpp, duty_cycle, clock_time, animal_id, rf_gain_dB,
                 photometry_ch1, photometry_ch2, file_names, file_index):
        # usp object for containing experimental parameters and data

        self.t_pre = t_pre  # time (s) before stimulus for each event
        self.t_on = t_on  # time (s) stimulus is on for each event
        self.t_post = t_post  # time (s) after stimulus in each event
        self.fs = fs  # sampling frequency for photometry in experiment
        self.f_c = f_c  # carrier frequency in experiment
        self.f_mod = f_mod  # freq modulation list; each item in array is f_mod for event
        self.vpp = vpp  # voltage list; each item in array is voltage for event
        self.duty_cycle = duty_cycle  # voltage list; each item in array is voltage for event
        self.clock_time = clock_time  # Clock date-time at beginning of experiment run
        self.animal_id = animal_id    # id of animal used
        self.rf_gain_db = rf_gain_dB  # gain from radio frequency amplifier in decibals
        self.file_names = file_names  # array of file names constituting usp object
        self.file_index = file_index  # array of indices for file_names association; useful with merged files
        self.vpp_filtered = []  # vpp event list filtered to include events containing specified vpp/f_mod
        self.f_mod_filtered = []  # f_mod event list filtered
        self.file_index_filtered = []  # file_index event list filtered
        self.subtracted_delta_f = []  # vector for channel subtraction data; ie isobestic signal subtraction
        self.fft_vals = []  # vector of fft values for each event of a given channel, can be overwritten when called
        self.fft_freqs = []   # vector of fft frequencies for each event of a given channel, can be overwritten
        self.phase_angles = []  # vector of phase angles at given point for each trial
        self.phase_degrees = []  # vector of phase degrees at given point for each trial
        self.photometry_ch = [photometry_ch1, photometry_ch2]  # 2D array with photometry channel 1 data for each event
        self.empty_chan_list = [[] for i in
                                range(len(self.photometry_ch))]  # Empty 2D array where length = # of channels
        self.delta_f = copy.deepcopy(self.empty_chan_list)  # data output for transforming data into delta values
        self.delta_f_filtered = copy.deepcopy(
            self.empty_chan_list)  # data output for transforming data into filtered delta values
        self.photometry_filtered = copy.deepcopy(
            self.empty_chan_list)  # photometry_ch1 filtered to include events containing specified vpp/f_mod
        self.slice_means = copy.deepcopy(self.empty_chan_list)  # get mean of time slice in raw photometry_ch1
        self.delta_f_slice_means = copy.deepcopy(self.empty_chan_list)  # get mean of time slice in delta_f
        self.rs = []

    def get_delta_f(self, baseline_start, baseline_end):
        # class function for converting photometry reads into delta f values using input baseline parameters
        # baseline_start - (s) starting time of baseline period
        # baseline_end   - (s) end time of baseline period

        if baseline_start < 0 or baseline_end > self.t_post + self.t_on + self.t_pre:
            print('Error: baseline values out of range. Entire trace may be used as baseline. ')

        for i in range(len(self.delta_f)):
            self.delta_f[i] = calculate_delta_f(self.photometry_ch[i], baseline_start, baseline_end, self.fs)

    def normalize_delta_f(self, channel=0, norm_base='abs_max'):
        # class function for converting photometry reads into delta f values using input baseline parameters
        # baseline_start - (s) starting time of baseline period
        # baseline_end   - (s) end time of baseline period

        if norm_base =='abs_max':
            # get max absolute delta value
            norm_val = get_abs_max_2d(self.delta_f[channel])

        if norm_base =='mean':
            # get mean delta value of all trial data
            norm_val = get_abs_mean_2d(self.delta_f[channel])  #this only works on type <class 'numpy.ndarray'>

        for i in range(len(self.delta_f[channel])):
            for val in range(len(self.delta_f[channel][i])):
                self.delta_f[channel][i][val] = self.delta_f[channel][i][val]/ norm_val

    def acquire_usp_data(self, usp_class_instance):
        # usp class method for merging data between two instances of usp class
        # requires that t_pre, t_on, t_post, f_s, f_c match between instances.
        # Instance name is inherited from the instance calling the method rather than the acquired instance
        # File indices attribute (file_index) of acquired file is incremented up by the number of files
        # contained in the acquiring file. File names stored in file_names attributes

        if usp_class_instance.t_pre == self.t_pre and usp_class_instance.t_on == self.t_on and \
                usp_class_instance.t_post == self.t_post and usp_class_instance.fs == self.fs:

            for i in range(len(self.photometry_ch)):
                self.photometry_ch[i] = np.concatenate((self.photometry_ch[i], usp_class_instance.photometry_ch[i]),
                                                       axis=0)
            self.f_mod = np.concatenate((self.f_mod, usp_class_instance.f_mod), axis=0)
            self.vpp = np.concatenate((self.vpp, usp_class_instance.vpp), axis=0)
            self.file_index = np.concatenate((self.file_index, usp_class_instance.file_index + len(self.file_names)),
                                             axis=0)
            self.file_names = self.file_names + usp_class_instance.file_names

    def filter_data(self, filter_vals=None, must_contain_all=False):
        # filter_data iterates through the vpp and f_mod array retrieving only events which
        # contain at least one item from filter_vals, returning a filtered copy of photometry
        # data, vpp, and f_mod list

        for i in range(len(self.photometry_ch)):
            self.photometry_filtered[i] = np.zeros((len(self.photometry_ch[i]), len(self.photometry_ch[i][0])))
        self.vpp_filtered = np.zeros(len(self.vpp))
        self.f_mod_filtered = np.zeros(len(self.f_mod))
        self.file_index_filtered = np.zeros(len(self.file_index))

        if must_contain_all == False:
            for i in range(len(self.vpp)):
                if self.vpp[i] in filter_vals or self.f_mod[i] in filter_vals:
                    for ch in range(len(self.photometry_filtered)):
                        self.photometry_filtered[ch][i] = self.photometry_ch[ch][i]
                    self.vpp_filtered[i] = self.vpp[i]
                    self.f_mod_filtered[i] = self.f_mod[i]
                    self.file_index_filtered = self.file_index[i]

        else:
            for i in range(len(self.vpp)):
                if self.vpp[i] in filter_vals and self.f_mod[i] in filter_vals:
                    for ch in range(len(self.photometry_filtered)):
                        self.photometry_filtered[ch][i] = self.photometry_ch[ch][i]
                    self.vpp_filtered[i] = self.vpp[i]
                    self.f_mod_filtered[i] = self.f_mod[i]
                    self.file_index_filtered = self.file_index[i]

        for i in range(len(self.photometry_filtered)):
            self.photometry_filtered[i] = self.photometry_filtered[i][~np.all(self.photometry_filtered[i] == 0, axis=1)]

        self.vpp = self.vpp_filtered[self.vpp_filtered != 0]
        self.f_mod = self.f_mod_filtered[self.f_mod_filtered != 0]
        self.file_index = self.file_index_filtered[self.file_index_filtered != 0]
        self.photometry_ch = self.photometry_filtered

        if len(self.vpp_filtered) == 0:
            print('Warning: given filter values do not exist in data')

    def subtract_delta_f(self, ch_subtracted_from, ch_subtracted):

        if self.delta_f != self.empty_chan_list:
            for i in range(len(self.photometry_ch[ch_subtracted_from])):
                self.subtracted_delta_f.append(
                    np.subtract(self.delta_f[ch_subtracted_from][i], self.delta_f[ch_subtracted][i]))

            self.delta_f[ch_subtracted_from] = self.subtracted_delta_f

        else:
            print('Delta calculation required!')

    def remove_decay(self, channel, padding=False, decay_type='exponential'):
        # this function will linearize all events to fit a sort of global exponential decay
        # since single events can be very short in which case a fit would not be obvious

        trace_set = get_whole_trace(self)
        pad = min(get_whole_trace(self)[channel])
        if decay_type == 'exponential':
            try:
                subtract_y = remove_exp_decay(trace_set[channel])
            except:
                print('Warning: Exponential decay fit failed. Using linear decay fit instead.')
                decay_type = 'linear'

        if decay_type == 'linear':
            subtract_y = remove_lin_decay(trace_set[channel])

        linearized_index = 0
        for trial in range(len(self.photometry_ch[channel])):
            for read in range(len(self.photometry_ch[channel][trial])):
                self.photometry_ch[channel][trial][read] = self.photometry_ch[channel][trial][read] - subtract_y[linearized_index]
                linearized_index = linearized_index + 1

        if padding == True:
            for trial in range(len(self.photometry_ch[channel])):
                for read in range(len(self.photometry_ch[channel][trial])):
                    self.photometry_ch[channel][trial][read] = self.photometry_ch[channel][trial][read] + pad

    def bin_data(self, bin_size):
        # bin photometry data by bin size in reads, edits fs value

        self.fs = self.fs/bin_size
        for i in range(len(self.photometry_ch)):

            temp_photometry_ch = np.zeros((self.photometry_ch[i].shape[0], int(self.photometry_ch[i].shape[1]/bin_size)))
            for n in range(len(self.photometry_ch[i])):
                temp_photometry_ch[n] = bin_vector(self.photometry_ch[i][n], bin_size, numpy_array=True)
            self.photometry_ch[i] = temp_photometry_ch

    def get_fft(self, channel, window_start, window_end):
        # window start and end in seconds

        for i in range(len(self.delta_f[channel])):
            fft_freq, fft_values = fft_hann(self.delta_f[channel][i][int(window_start * self.fs):int(window_end * self.fs)], self.fs)
            self.fft_vals.append(fft_values)
            self.fft_freqs.append(fft_freq)
            #plt.plot(fft_freq, fft_values, color=[.7, 0, 0.5], alpha=0.5, linewidth=2)
            #plt.ylabel('Amplitude')
            #plt.xlabel('Frequency [Hz]')
            #plt.xlim(0,200)
            #plt.show()

        self.phase_angles, self.phase_degrees = get_stimulus_phase(self, window_start, window_end, plot_events=False)

    def get_slice_means(self, slice_start, slice_end):

        for i in range(len(self.photometry_ch)):
            self.slice_means[i] = np.zeros(len(self.photometry_ch[i]))
            self.delta_f_slice_means[i] = np.zeros(len(self.delta_f[i]))

        for i in range(len(self.photometry_ch[0])):
            for ch in range(len(self.photometry_ch)):
                self.slice_means[ch][i] = np.mean(
                    self.photometry_ch[ch][i][int(slice_start * self.fs):int(slice_end * self.fs)])
                self.delta_f_slice_means[ch][i] = np.mean(
                    self.delta_f[ch][i][int(slice_start * self.fs):int(slice_end * self.fs)])

    def write_file(self, file_name, swap_channels=False):
        # must be rewritten for different channel counts

        if swap_channels == False:
            dict_out = {'t_pre': self.t_pre, 't_on': self.t_on, 't_post': self.t_post, 'fs': self.fs, 'f_c': self.f_c,
                        'f_mod': self.f_mod, 'vpp': self.vpp, 'duty_cycle': self.duty_cycle, 'clock_time': self.clock_time,
                        'photometry_ch1': self.photometry_ch[0], 'photometry_ch2': self.photometry_ch[1],
                        'file_names': self.file_names, 'file_index': self.file_index, 'animal_id': self.animal_id,
                        'rf_gain_db': self.rf_gain_db}

        if swap_channels == True:
            dict_out = {'t_pre': self.t_pre, 't_on': self.t_on, 't_post': self.t_post, 'fs': self.fs, 'f_c': self.f_c,
                        'f_mod': self.f_mod, 'vpp': self.vpp, 'duty_cycle': self.duty_cycle, 'clock_time': self.clock_time,
                        'photometry_ch1': self.photometry_ch[1], 'photometry_ch2': self.photometry_ch[0],
                        'file_names': self.file_names, 'file_index': self.file_index, 'animal_id': self.animal_id,
                        'rf_gain_db': self.rf_gain_db}

        np.save(file_name, dict_out)

    def remove_events(self, event_ch_1, event_ch_2, criteria=0):
        # Should be done before any calculations on data are made
        if self.empty_chan_list != self.delta_f:
            print('error: must be performed prior to data analysis')
            return

        delete_index = []
        mean_list = []
        for i in range(len(self.photometry_ch[event_ch_1])):
            # if some criteria are true, add index to be deleted
            base_mean = np.mean(self.photometry_ch[event_ch_1][i][0:100])
            mean_list.append(base_mean)
            if i > 30:
                delete_index.append(i)
        #usp_graphics.histogram(mean_list)

        print("Removing " + str(len(delete_index) / len(self.photometry_ch[event_ch_1]) * 100) + "% of events")
        for n in range(len(self.photometry_ch)):
            self.photometry_ch[n] = np.delete(self.photometry_ch[n], delete_index, axis=0)
        self.clock_time = np.delete(self.clock_time, delete_index, axis=0)
        self.vpp = np.delete(self.vpp, delete_index, axis=0)
        self.f_mod = np.delete(self.f_mod, delete_index, axis=0)

    def remove_events_by_phase(self, event_ch_1, event_ch_2, criteria=0):
        # must be done following usp.fft analysis

        delete_index = []

        for i in range(len(self.photometry_ch[event_ch_1])):
            # if some criteria are true, add index to be deleted

            if self.phase_degrees[i] < 180 or self.phase_degrees[i] > 270:
                # if base_mean > 5.1:
                delete_index.append(i)

        print("Removing " + str(len(delete_index) / len(self.photometry_ch[event_ch_1]) * 100) + "% of events")
        for n in range(len(self.photometry_ch)):
            self.photometry_ch[n] = np.delete(self.photometry_ch[n], delete_index, axis=0)
            self.delta_f[n] = np.delete(self.delta_f[n], delete_index, axis=0)
        self.clock_time = np.delete(self.clock_time, delete_index, axis=0)
        self.vpp = np.delete(self.vpp, delete_index, axis=0)
        self.f_mod = np.delete(self.f_mod, delete_index, axis=0)

    def subtract_out(self, sub_value):

        for i in range(len(self.photometry_ch)):
            for trial in range(len(self.photometry_ch[i])):
                for read in range(len(self.photometry_ch[i][trial])):
                    self.photometry_ch[i][trial][read] = self.photometry_ch[i][trial][read] - sub_value

    def smooth_signal(self, channel, window=20):
        # must be rewritten for different channel counts
        for i in range(len(self.photometry_ch[channel])):
            self.photometry_ch[channel][i] = moving_average(self.photometry_ch[channel][i], window)

    def adjust_gain_delta_f(self, channel, gain=1):

        for i in range(len(self.delta_f[channel])):
            self.delta_f[channel][i] = self.delta_f[channel][i] * gain
