from tkinter import *
from tkinter import filedialog
import os
import mne
import numpy as np
import matplotlib.pyplot as plt


def browseFiles():
    cdpath = os.getcwd()
    filename = filedialog.askopenfilename(initialdir=cdpath, title="Select a raw file",
                                          filetypes=[("MNE data files", "*.*")])
    return filename


def _dolowpass(rawdata, sfreq_aim, chnoms_all, chnoms):
    """Here we low-pass the data at 0.33 of the desired sample rate.
       This is a first phase of downsampling, which aims to avoid
       the reduction of temporal precision of events often introduced
       with resampling.
       Decimation will be carried out on the segmented data.
       Note: For this approach to function,the original sampling frequency must be
       an integer multiple of the new sampling frequency.
       """
    sfreq_curr = rawdata.info["sfreq"]
    decim = np.round(sfreq_curr / sfreq_aim).astype(int)
    sfreq_new = sfreq_curr / decim
    freq_lowpass = sfreq_new / 3
    print('The low-pass frequency is {} Hz'.format(freq_lowpass))

    chindx = [chnoms_all.index(ic) for ic in chnoms]  # Find the indices of the 'eeg' channels.
    trans_bandwidth = min(max(freq_lowpass * 0.25, 2.),
                          sfreq_curr / 2. - freq_lowpass)  # Find the transition bandwidth (Hz) for upper cutoff.
    rawfilt_lp = rawdata.copy().filter(l_freq=None, h_freq=freq_lowpass, picks=chindx, filter_length='auto',
                                       h_trans_bandwidth=trans_bandwidth, n_jobs=1, fir_design='firwin')
    return rawfilt_lp


def _dohighpass(rawdata, hpfreq, chnoms_all, chnoms):
    """Carry out high-pass filtering on the low-pass filtered data."""

    chindx = [chnoms_all.index(ic) for ic in chnoms]
    trans_bandwidth = min(max(hpfreq * 0.25, 2), hpfreq)

    rawfilt_hp = rawdata.copy().filter(l_freq=hpfreq, h_freq=None, picks=chindx, filter_length ='auto',
                                       l_trans_bandwidth=trans_bandwidth, n_jobs=1, fir_design ='firwin')
    return rawfilt_hp

def plotSpectrum(rawdata):
    """Function to plot the spectrum of selected electrodes.
       This can be useful for detecting noisy electrodes.
       """
    psd_out, freq_out = mne.time_frequency.psd_welch(rawdata, fmin=0.4, fmax=80, tmin=100, tmax=400, picks=[])


"""****************************** LOAD IN CONTINUOUS DATA IN .FIF FORMAT################################"""
filename = browseFiles()
rawIn = mne.io.read_raw_fif(filename, allow_maxshield=True, preload=True, verbose=None)

# Extract basic information from the rawIn.info
# Need to separate the channel names from the trigger channel names, which begin with "D"
ch_names_all = rawIn.ch_names
ch_names = [idx for idx in ch_names_all if idx.startswith("E")]
misc_names = [idx1 for idx1 in ch_names_all if idx1.startswith("D")]

rawIn.set_channel_types(dict.fromkeys(ch_names, 'eeg'))  # Set the ch_names to 'eeg' type.
rawIn.set_channel_types(dict.fromkeys(misc_names, 'misc'))  # Set the trig channel type to 'misc'
print('The eeg scalp channels are as follows: \n')
print(rawIn.copy().pick_types(meg=False, eeg=True).ch_names)
print('The trig channels are as follows: \n')
print(rawIn.copy().pick_types(meg=False, misc=True).ch_names)
hpfilt = rawIn.info['highpass']
print("High-pass filter cutoff: " + str(hpfilt) + "Hz")  # Print high-frequency cut-off frequency to screen.
lowfilt = rawIn.info['lowpass']
print("Low-pass filter cutoff: " + str(lowfilt) + "Hz")  # Print low-frequency cut-off frequency to screen.

# Exract the time vector in seconds
time_secs = rawIn.times

print('Sampling frequency of current dataset: ', rawIn.info['sfreq'], 'Hz', '\n')
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

Events = mne.find_events(rawIn)  # Extract the events as an ndarray
raw_eeg = rawIn.copy().pick_types(meg=False, eeg=True)  # Only show the EEG channels.
mne.viz.plot_raw(raw_eeg, events=Events, duration=10.0, n_channels=10, title='Raw EEG Data', event_color='red',
                 remove_dc=True)

### ************* Set the electrode montage for the current data ************************
"""Here using GSN-HydroCel-257 montage. 
   Visualizing the montage in 2D.
   Assign montage to raw object (rawIn).
   """
montage_all = mne.channels.get_builtin_montages()
montindx = montage_all.index('GSN-HydroCel-257')
montage = mne.channels.make_standard_montage(montage_all[montindx])
montage.rename_channels({'Cz': 'E257'})  # Changing channel 'Cz' in montage to 'E257.
mne.viz.plot_montage(montage)  # Visualize montage.
rawIn.set_montage(montage)  # Assign the montage to the rawIn object.

##*********************** Carry out low-pass and high-pass filtering ******************
new_sfreq = 250                                                     # Define the sampling frequency (Hz) desired after downsampling.
rawfilt_LP = _dolowpass(rawIn, new_sfreq, ch_names_all, ch_names)   # Call of function to carry out low-pass filtering.

hplim = .1                                                          # High-pass filter cutoff, .1Hz
rawfilt_HP = _dohighpass(rawfilt_LP, hplim, ch_names_all, ch_names) # Call of function to carry out high-pass filtering.

### ********************* Detecting Bad Electrodes *************************
# Another way of calculating and plotting the PSD in MNE.
