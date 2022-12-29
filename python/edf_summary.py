'''
This code read the EDF file and creates summary details
'''

import mne
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mne
import glob
from numpy import loadtxt

fname = '//Users/sq566/Downloads/nsrr_studies/apples/subjects.txt'
subjects = loadtxt(fname, dtype='str').tolist()

subject_name = []
frequency_resolution = []
last_time_sample = []
sampling_frequency = []
raw_data_shape = []
sampling_rate_same_across_ch = []
ch_names_list = []

for sub in subjects:
    edf_file = f'/Users/sq566/Downloads/nsrr_studies/apples/data/edf/{sub}.edf'
    raw = mne.io.read_raw_edf(edf_file, preload=False,verbose=False)
    data = raw._raw_extras[0]
    ch_names = data['ch_names']
    meas_date = data['meas_date']
    n_samps = data['n_samps']
    offset = data['offsets']
    times = raw.times
    
    # Check if all the channels have the same sampling frequency rate
    result = all(element == n_samps[0] for element in n_samps) 
    
    srate = raw.info['sfreq']
    n_time_samps = raw.n_times
    time_secs = raw.times
    ch_names = raw.ch_names
    n_chan = len(ch_names)
    freq_res =  srate/n_time_samps

    subject_name.append(sub)
    frequency_resolution.append(freq_res)
    last_time_sample.append(round(time_secs[-1]))
    sampling_frequency.append(srate)
    raw_data_shape.append(f'({n_chan} x {n_time_samps})')
    sampling_rate_same_across_ch.append(result)
    ch_names_list.append(ch_names)
    

info_dict = {'subject': subject_name, 'last_time_sample': last_time_sample, 
              'sampling_frequency':sampling_frequency, 'raw_data_shape':raw_data_shape, 'sampling_rate_same_across_ch':sampling_rate_same_across_ch,'ch_names':ch_names_list}

df = pd.DataFrame(info_dict)
output_fname = '/Users/sq566/Downloads/nsrr_studies/apples/info_summary.xlsx'
df.to_excel(output_fname, index=False)
