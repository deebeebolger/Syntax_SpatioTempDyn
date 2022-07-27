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
rawIn = mne.io.read_raw_fif(filename, allow_maxshield=True, preload=False, verbose=None)

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


# Exract the time vector in seconds
time_secs     = rawIn.times

print('Sampling frequency of current dataset: ',rawIn.info['sfreq'],'Hz','\n')
print('The current dataset has {} time samples and {} channels \n'.format(rawIn.n_times, len(ch_names)))
print('The duration of the current dataset is: {}seconds'.format(time_secs[-1]))

# ************ Simple Visualization of all Channels **********************************
""" Note that in the following we are only visualizing the EEG channels
    The time interval (in seconds) and the number of electrodes presented in the main window is defined
    when calling the visualisation method.
    The events are indicated on the plot (blue lines).
    The DC offset is applied for visualisation purposes because is not yet high-pass filtered.
    You change the scale of visualisation using the '+' and '-' keys.
    Time intervals with data missing are/should be indicated on the time-bar at the bottom of the figure.
    Here you can begin to identify bad/noisy electrodes."""

Events = mne.find_events(rawIn)
raw_eeg = rawIn.copy().pick_types(meg=False, eeg=True)                      # Only show the EEG channels.
mne.viz.plot_raw(raw_eeg, events=Events, duration=10.0, n_channels=10,title='Raw EEG Data', remove_dc=True)




