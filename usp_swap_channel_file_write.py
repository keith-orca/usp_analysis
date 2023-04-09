import usp_load_data as lmd
from usp_graphics import *
from usp_functions import *

file_name = 'bns6_20hz_40s_intensity'
usp = lmd.load_npy_file('C:\\Users\\15163\\Desktop\\stanford_software\\usp_analysis_code\\usp_data\\' + file_name + '.npy')

# Write the usp object to a numpy dictionary file; this is useful if files have been merged
usp.write_file(file_name, swap_channels=True)
