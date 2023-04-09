from scipy.io import loadmat
from usp_us_photometry_class import *
import os


def load_npy_file(npy_file_name):
    # Loads npy dictionary file containing us photometry class data
    # This file is typically generated when two mat struct files are
    # merged and so file_names will reflect constituents

    data = np.load(npy_file_name, allow_pickle=True)
    data = data[()]  # gets dict from numpy.ndarray object
    # Below are of the class int
    t_pre = data.get('t_pre')
    t_on = data.get('t_on')
    t_post = data.get('t_post')
    fs = data.get('fs')
    f_c = data.get('f_c')
    # Below are of the class numpy.ndarray
    f_mod = data.get('f_mod')
    vpp = data.get('vpp')
    duty_cycle = data.get('duty_cycle')
    clock_time = data.get('clock_time')
    photometry_ch1 = data.get('photometry_ch1')
    photometry_ch2 = data.get('photometry_ch2')
    file_names = data.get('file_names')
    file_index = data.get('file_index')
    animal_id = data.get('animal_id')
    rf_gain_db = data.get('rf_gain_dB')
    us_photometry_instance = USPhotometry(t_pre, t_on, t_post, fs, f_c, f_mod, vpp, duty_cycle, clock_time, animal_id,
                                          rf_gain_db, photometry_ch1, photometry_ch2, file_names, file_index)

    return us_photometry_instance


def load_folder(directory_name):
    # scan folder and merge all usp files

    files = []
    file_id = []
    usp_set = []

    for file_name in os.listdir(directory_name):
        if file_name.find('.mat') != -1:
            print('reading ' + file_name + '...')
            usp_set.append(load_mat_file(os.path.join(directory_name, file_name)))

        if file_name.find('.npy') != -1:
            print('reading ' + file_name + '...')
            usp_set.append(load_npy_file(os.path.join(directory_name, file_name)))

    for i in range(1, len(usp_set)):
        usp_set[0].acquire_usp_data(usp_set[i])

    return usp_set[0]


def load_mat_file(mat_file_name):
    # No longer in use...
    # function for importing .mat file as dictionary and loading parameters into usp class object

    data = loadmat(mat_file_name, squeeze_me=True, struct_as_record=False)

    # parameters of the class int
    t_pre = data['results'].t_pre
    t_on = data['results'].t_on
    t_post = data['results'].t_post
    fs = data['results'].fs
    f_c = data['results'].f_c

    # parameters are of the class numpy.ndarray
    f_mod = data['results'].f_mod
    vpp = data['results'].Vpp
    duty_cycle = None
    clock_time = data['results'].clock_time
    photometry_ch1 = data['results'].photometry_CH1.transpose()
    photometry_ch2 = data['results'].photometry_CH2.transpose()
    file_names = [strip_file_path(mat_file_name)]
    file_index = np.zeros(len(vpp))
    animal_id = ''
    rf_gain_dB = ''

    # instantiate the us photometry instance from file parameters
    us_photometry_instance = USPhotometry(t_pre, t_on, t_post, fs, f_c, f_mod, vpp, duty_cycle, clock_time, animal_id,
                                          rf_gain_dB, photometry_ch1, photometry_ch2, file_names, file_index)

    return us_photometry_instance