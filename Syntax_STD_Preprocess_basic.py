from tkinter import *
from tkinter import filedialog
import os
import mne
import numpy as np
import matplotlib.pyplot as plt

def browseFiles():
    cdpath = os.getcwd()
    filename = filedialog.askopenfilename(initialdir = cdpath, title="Select a raw file", filetypes=[("MNE data files", "*.*")])
    return filename

"""****************************** LOAD IN CONTINUOUS DATA IN .FIF FORMAT################################"""
filename = browseFiles()
rawIn = mne.io.read_raw_fif(filename, allow_maxshield=True, preload=True, verbose=None)

# Extract basic information from the rawIn.info
# Need to separate the channel names from the trigger channel names, which begin with "D"
ch_names_all  = rawIn.ch_names
ch_names      = [idx for idx in ch_names_all if idx.startswith("E")]
misc_names    = [idx1 for idx1 in ch_names_all if idx1.startswith("D")]

rawIn.set_channel_types(dict.fromkeys(ch_names,'eeg'))       # Set the ch_names to 'eeg' type.
rawIn.set_channel_types(dict.fromkeys(misc_names,'misc'))    # Set the trig channel type to 'misc'
print('The eeg scalp channels are as follows: \n')
print(rawIn.copy().pick_types(meg=False, eeg=True).ch_names)
print('The trig channels are as follows: \n')
print(rawIn.copy().pick_types(meg=False, misc=True).ch_names)
hpfilt = rawIn.info['highpass']
print("High-pass filter cutoff: "+str(hpfilt)+"Hz")      #Print high-frequency cut-off frequency to screen.
lowfilt = rawIn.info['lowpass']
print("Low-pass filter cutoff: "+str(lowfilt)+"Hz")      #Print low-frequency cut-off frequency to screen.

# Exract the time vector in seconds
time_secs     = rawIn.times

print('Sampling frequency of current dataset: ',rawIn.info['sfreq'],'Hz','\n')
print('The current dataset has {} time samples and {} channels \n'.format(rawIn.n_times, len(ch_names)))
print('The duration of the current dataset is: {}seconds'.format(time_secs[-1]))

# ************ Simple Visualization of all Channels **********************************
""" Note that in the following we are only visualizing the EEG channels
    The time interval (in seconds) and the number of electrodes presented in the main window is defined
    when calling the visualisation method.
    The events are indicated on the plot (red lines).
    The DC offset is applied for visualisation purposes because is not yet high-pass filtered.
    You change the scale of visualisation using the '+' and '-' keys on your keyboard.
    You can scroll through the data in time using the right/left arrows on your keyboard.
    Time intervals with data missing are/should be indicated on the time-bar at the bottom of the figure.
    Here you can begin to identify bad/noisy electrodes."""

Events = mne.find_events(rawIn)                                       # Extract the events as an ndarray
raw_eeg = rawIn.copy().pick_types(meg=False, eeg=True)                # Only show the EEG channels.
mne.viz.plot_raw(raw_eeg, events=Events, duration=10.0, n_channels=10,title='Raw EEG Data', event_color='red', remove_dc=True)

### ************* Set the electrode montage for the current data ************************
"""Here using GSN-HydroCel-257 montage. 
   Visualizing the montage in 2D.
   """
montage_all = mne.channels.get_builtin_montages()
montindx = montage_all.index('GSN-HydroCel-257')
montage = mne.channels.make_standard_montage(montage_all[montindx])
montage.rename_channels({'Cz':'E257'})
mne.viz.plot_montage(montage)
rawIn.set_montage(montage)  # Assign the montage to the rawIn object.

#%% Downsample the continuous data from 500Hz to 250Hz ***************************
newfreq = 250
rawrs = rawIn.copy().resample(sfreq=newfreq)
print('The new sampling frequency is: ',rawrs.info['sfreq'],'Hz')

#%% Apply a FIR bandpass filter which will conserve activity in the 0.1Hz to 80Hz frequency band and reject all others.
"""As filtering transforms the raw object, we make a copy of it and apply the filter to this copy; 
   rawfilt is now our filtered object.
   We are filtering the downsampled data (rawrs).
   The data is band pass filtered between 0.1Hz and 80Hz.
   We can see if we can push the lower limit down to 0.05Hz without leaving too 
   much slow drift.
   Note that the raw object data has to be preloaded upon import to be filtered.
   """
rawfilt = rawrs.copy().filter(0.1, 80, n_jobs=2, fir_design='firwin')   # need to look up the filter parameters again!!

### ********************* Detecting Bad Electrodes *************************
# Another way of calculating and plotting the PSD in MNE.

psd_out, freq_out = mne.time_frequency.psd_welch(rawfilt, fmin=0.4, fmax=80, tmin=100, tmax=400, picks=[])




