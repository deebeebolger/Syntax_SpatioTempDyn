from tkinter import *
from tkinter import filedialog
import os
import mne
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from meegkit import dss
from meegkit.utils import create_line_data, unfold

def browseFiles():
    cdpath = os.getcwd()
    filename = filedialog.askopenfilename(initialdir=cdpath, title="Select a raw file",
                                          filetypes=[("MNE data files", "*.*")])
    return filename

def viz_allchans(rawdata, vizdur, channum_viz):
    """ Simple Visualization of all Channels stacked
        Note that in the following we are only visualizing the EEG channels
        The time interval (in seconds) and the number of electrodes presented in the main window is defined
        when calling the visualisation method.
        The events are indicated on the plot (red lines).
        The DC offset is applied for visualisation purposes because is not yet high-pass filtered.
        You change the scale of visualisation using the '+' and '-' keys on your keyboard.
        You can scroll through the data in time using the right/left arrows on your keyboard.
        Time intervals with data missing are/should be indicated on the time-bar at the bottom of the figure.
        Here you can begin to identify bad/noisy electrodes by simply selecting them on the figure (just clicking on them).
        When you close the visualisation window the channels you marked as bad will be displayed in the console and will appear in
        rawdata.info['bads]
        """

    Events = mne.find_events(rawdata)  # Extract the events as an ndarray
    raw_eeg = rawdata.copy().pick_types(meg=False, eeg=True)  # Only show the EEG channels.
    figviz = mne.viz.plot_raw(raw_eeg, events=Events, duration=vizdur, n_channels=channum_viz, title='Raw EEG Data', event_color='red',
                     remove_dc=True, block=False, show=True)

    bads = raw_eeg.info['bads']
    return bads

def _dohighpass(rawdata, hpfreq, chnoms_all, chnoms):
    """Carry out high-pass filtering on the low-pass filtered data."""

    chindx = [chnoms_all.index(ic) for ic in chnoms]
    trans_bandwidth = min(max(hpfreq * 0.25, 2), hpfreq)
    filter_params = mne.filter.create_filter(rawdata.get_data(), rawdata.info["sfreq"], l_freq=hpfreq, h_freq=None,
                                             l_trans_bandwidth=trans_bandwidth,  fir_design='firwin2',
                                             verbose=True)  # Acquire filter parameters
    hpfilt_fig = mne.viz.plot_filter(filter_params, rawdata.info["sfreq"])
    rawfilt_hp = rawdata.copy().filter(l_freq=hpfreq, h_freq=None, picks=chindx, filter_length ='auto',
                                       l_trans_bandwidth=trans_bandwidth, n_jobs=1, fir_design ='firwin2')
    return rawfilt_hp, hpfilt_fig, filter_params

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
    freq_lowpass = sfreq_new / 3        # Lowpass filtering at the 0.33 of the desired sample rate.
    print('The low-pass frequency is {} Hz'.format(freq_lowpass))

    chindx          = [chnoms_all.index(ic) for ic in chnoms]                                                       # Find the indices of the 'eeg' channels.
    trans_bandwidth = min(max(freq_lowpass * 0.25, 2.), sfreq_curr / 2. - freq_lowpass)                             # Find the transition bandwidth (Hz) for upper cutoff.
    filter_params   = mne.filter.create_filter(rawdata.get_data(), sfreq_curr, h_freq=freq_lowpass,
                                               h_trans_bandwidth=trans_bandwidth, l_freq=None,
                                               fir_design='firwin2', verbose=True)  # Acquire filter parameters
    lpfilt_fig      = mne.viz.plot_filter(filter_params, sfreq_curr)
    rawfilt_lp      = rawdata.copy().filter(l_freq=None, h_freq=freq_lowpass, picks=chindx, filter_length='auto',
                                       h_trans_bandwidth=trans_bandwidth, n_jobs=1, fir_design='firwin2')
    return rawfilt_lp, lpfilt_fig, filter_params



def findNextPowerOf2(L):
    # decrement `n` (to handle cases when `n` itself
    # is a power of 2)
    L = L - 1

    # do till only one bit is left
    while L & L - 1:
        L = L & L - 1  # unset rightmost bit

    # `n` is now a power of two (less than `n`)

    # return next power of 2
    return L << 1

def plotSpectrum(rawdata, powtwoL, chindx):
    """Function to plot the spectrum of selected electrodes.
       This can be useful for detecting noisy electrodes.
    """
    nfft_len = powtwoL
    L = rawdata.last_samp
    psd_out, freq_out = mne.time_frequency.psd_welch(rawdata, fmin=0.4, fmax=60, tmin=0, tmax=rawdata.times.max(), picks=chindx, n_fft=256, reject_by_annotation=TRUE)
    np.seterr(divide='ignore')
    log_psd = 10*np.log10(psd_out)
    np.seterr(divide='warn')
    _, ax = plt.subplots()
    ax.plot(freq_out, log_psd.T)
    ax.set(title='Welch PSD of band-pass filtered data', xlabel='Frequency (Hz)', ylabel='Power Spectral Density (dB)')
    plt.show()

def create_savepath(fstring, sfix):
    fs = fstring.split('/')
    datacurr = fs[-1]
    data2save = datacurr.split('.')
    rawout_name = data2save[0] + '_' + sfix + '.fif'
    fs1 = fs[0:-1]
    savepath = '/'.join(fs1)
    savepath = savepath + '/'
    return savepath, rawout_name

"""****************************** LOAD IN CONTINUOUS DATA IN .FIF FORMAT################################"""
study_name = 'DynSyn_EEG_PC'
study_folder = '_'.join(['Project', study_name])
study_root =  os.path.join("Users", "bolger", "Documents", "work", "Projects", study_folder)
bids_root  =  os.path.join(study_root, study_name)    # BIDS root directory

sessions = ['1']
subjects = ['05']
ch_types = ['eeg']
data_type = ('eeg')

## ************* Define the directory in which to search for the raw.fif data *****************
subdir = '-'.join(['sub',subjects[0]])
sesdir = '-'.join(['ses', sessions[0]])
bids_suj = os.path.join(bids_root, subdir, sesdir, data_type)     # The directory from which to load the current subject dataset.
if os.path.exists(bids_suj)==FALSE:
    raise Exception("The necessary directory for the current subject does not exist")
else:
    print('Loading raw file from subject directory: {}\n'.format(bids_suj))

## *********** If the "derivatives" directory does not already exist, create one
deriv_root = os.path.join(bids_root, 'derivatives')
if os.path.exists(deriv_root)==FALSE:
    os.mkdir(deriv_root)
else:
    print('Directory: {} and path: {} already exists.\n'.format('derivatives', deriv_root))

# If a directory for the current subject does not already exist, create one.
deriv_suj_root = os.path.join(deriv_root, subdir, sesdir, data_type)
if os.path.exists(deriv_suj_root)==FALSE:
    os.mkdir(deriv_suj_root)
else:
    print('Directory: {} and path: {} already exists.\n'.format(subdir, deriv_suj_root))





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
print('High-pass filter cutoff: ', hpfilt, 'Hz')  # Print high-frequency cut-off frequency to screen.
lowfilt = rawIn.info['lowpass']
print('Low-pass filter cutoff: ', lowfilt,'Hz')  # Print low-frequency cut-off frequency to screen.

# Exract the time vector in seconds
time_secs = rawIn.times

print('Sampling frequency of current dataset: ', rawIn.info['sfreq'], 'Hz', '\n')
print('The current dataset has {} time samples and {} channels \n'.format(rawIn.n_times, len(ch_names)))
print('The duration of the current dataset is: {}seconds'.format(time_secs[-1]))

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
  # Call of function to carry out low-pass filtering.
new_sfreq = 250                                                     # Define the sampling frequency (Hz) desired after downsampling.
rawfilt_LP = _dolowpass(rawIn, new_sfreq, ch_names_all, ch_names)
hplim = .1                                                          # High-pass filter cutoff, .1Hz
rawfilt_HP = _dohighpass(rawfilt_LP, hplim, ch_names_all, ch_names) # Call of function to carry out high-pass filtering.

### ********************* Detecting Bad Electrodes *************************

lensig = rawfilt_HP.last_samp
round(len(rawfilt_HP)/10)
pow2 = findNextPowerOf2(int(lensig/2))
chanindx = [ch_names_all.index(ic) for ic in ch_names]
print('Signal length is {} and the nearest power of two is {}'.format(lensig, pow2))
plotSpectrum(rawfilt_HP, pow2, chanindx)
print('****plotted spectrum')

### ******* Apply the zpline function (from the meegkit package) to remove line noise *********####

## Call of function here in which the zapline function is run.
dss.dss_line(rawfilt_HP, 50, new_sfreq, nfft=500)

## ********************* Call of function viz_allchans() to plot all chans
wind_dur    = 10         # Duration of time interval (in seconds) presented in each window.
wind_nchans = 50         # Number of channels presented in each window.
badchans = viz_allchans(rawfilt_HP, wind_dur, wind_nchans)
print('The channels pre-selected as bad are: {}'.format(badchans))

suffix = 'bandpass'
path2save, fsave_name = create_savepath(filename,suffix)
rawfilt_HP.save(path2save + fsave_name)
