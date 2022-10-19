from tkinter import *
from tkinter import filedialog
import os
import mne
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import json
from autoreject import Ransac

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
    freq_lowpass = sfreq_new / 6        # Lowpass filtering at the 0.33 of the desired sample rate.
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


# def zapline_apply():
#     """ Apply the zapline function (from the meegkit package) to remove line noise
#         Call of function here in which the zapline function is run.
#     """
#     Events = mne.find_events(Rawfilt_LP)
#     mne.viz.plot_raw(RawIn, events=Events, duration=20, n_channels=20, title='Select channels with line noise',
#                      event_color='red',
#                      remove_dc=True, block=False, show=True)
#     badchans_line = Rawfilt_LP.info['bads']
#     chans_indx = mne.pick_channels(ch_names, include=badchans_line)  # Extract the indices of the channels of interest
#     Rdata1, times = Rawfilt_LP[chans_indx, 0:len(time_secs)]
#     Rdata1_tp = np.transpose(Rdata1)
#     Rawdss, artif = dss.dss_line_iter(Rdata1_tp, 50, sfreq, nfft=1000, show=True)
#     Rawdss_tp = np.transpose(Rawdss)
#     artif_tp = np.transpose(artif)
#
#     plt.plot(time_secs, Rdata1[10, :])
#     plt.plot(time_secs, Rawdss_tp[10, :])
#     plt.plot(time_secs, artif_tp[10, :])
#
#     for chindx, chs in enumerate(chans_indx):
#         Rawfilt_LP[chs, :] = np.repeat(0, len(time_secs))
#         Rawfilt_LP[chs, :] = Rawdss_tp[chindx, :]

"""****************************** LOAD IN CONTINUOUS DATA IN .FIF FORMAT################################"""
study_name = 'DynSyn_EEG_PC'
study_folder = '_'.join(['Project', study_name])
study_root =  os.path.join(" ", "Users", "bolger", "Documents", "work", "Projects", study_folder)
study_root = study_root[1:]
bids_root  =  os.path.join(study_root, study_name)    # BIDS root directory

sessions = ['1']
subjects = ['05']
ch_types = ['eeg']
data_type = ('eeg')

"""************* Define the directory in which to search for the raw.fif data *****************"""

subdir = '-'.join(['sub',subjects[0]])
sesdir = '-'.join(['ses', sessions[0]])
sujname = '_'.join([subdir, sesdir])
bids_suj = os.path.join(bids_root, subdir, sesdir, data_type)     # The directory from which to load the current subject dataset.
if os.path.exists(bids_suj)==FALSE:
    raise Exception("The necessary directory for the current subject does not exist")
else:
    print('Loading raw file from subject directory: {}\n'.format(bids_suj))

"""  If the "derivatives" directory does not already exist, create one
     The processed files are saved in the derivatives folder. 
     For each subject, a directory is created. 
 """
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

#%% ****************** LOAD IN THE CURRENT SUBJECT RAW FILE ************************** %%#
script_path_all = os.path.abspath('DynSynEEGPC_preprocessing.py')
script_exl      = script_path_all.split('/')
script_exl_path = '/'.join(script_exl[:-1])
htmlrep_path    = os.path.join(script_exl_path, 'MNE_Reports_html')
if os.path.exists(htmlrep_path)==FALSE:
    os.mkdir((htmlrep_path))
else:
    print('Directory: {} and path: {} already exists.\n'.format('MNE_Reports_html', htmlrep_path))

x = os.listdir(bids_suj+'/')
currf         = [ix for ix in x if ix.endswith('fif')]
fullpath_bids = os.path.join(bids_suj, currf[0])
RawIn         = mne.io.read_raw_fif(fullpath_bids, allow_maxshield=True, preload=True, verbose=None)
RawIn.pick_types(eeg=True, stim=True).load_data()

#%% ***** Set up an MNE Report (HTML file) ***********************
report_fname = '.'.join([currf[0][:-4], 'html'])
report_path  = os.path.join('MNE_Reports_html', report_fname)
data_report = mne.Report(title=currf[0], subject=subdir)
data_report.add_raw (raw=RawIn, title= "Raw data", psd=False)
data_report.save(fname = report_path, overwrite=True)

#%% ******************* Extract basic information from the rawIn.info
# Need to separate the channel names from the trigger channel names, which begin with "D"
ch_names_all = RawIn.ch_names
ch_names   = [idx for idx in ch_names_all if idx.startswith("E")]
misc_names = [idx1 for idx1 in ch_names_all if idx1.startswith("D")]

RawIn.drop_channels(ch_names_all[256])
RawIn.set_channel_types(dict.fromkeys(ch_names, 'eeg'))  # Set the ch_names to 'eeg' type.
RawIn.set_channel_types(dict.fromkeys(misc_names, 'misc'))  # Set the trig channel type to 'misc'
print('The eeg scalp channels are as follows: \n')
print(RawIn.copy().pick_types(meg=False, eeg=True).ch_names)
print('The trig channels are as follows: \n')
print(RawIn.copy().pick_types(meg=False, misc=True).ch_names)

sfreq = RawIn.info['sfreq']
time_secs = RawIn.times                           # Extract the time vector in seconds
hpfilt = RawIn.info['highpass']
print('High-pass filter cutoff: ', hpfilt, 'Hz')  # Print high-frequency cut-off frequency to screen.
lowfilt = RawIn.info['lowpass']
print('Low-pass filter cutoff: ', lowfilt,'Hz')   # Print low-frequency cut-off frequency to screen.
print('Sampling frequency of current dataset: ', RawIn.info['sfreq'], 'Hz', '\n')
print('The current dataset has {} time samples and {} channels \n'.format(RawIn.n_times, len(ch_names)))
print('The duration of the current dataset is: {}seconds'.format(time_secs[-1]))


### ************* Set the electrode montage for the current data ************************
"""Here using GSN-HydroCel-256 montage. 
   Visualizing the montage in 2D.
   Assign montage to raw object (rawIn).
"""
montage_all = mne.channels.get_builtin_montages()
montindx = montage_all.index('GSN-HydroCel-256')
montage = mne.channels.make_standard_montage(montage_all[montindx])
montage_fig = mne.viz.plot_montage(montage, show=False)  # Visualize montage.
RawIn.set_montage(montage)  # Assign the montage to the rawIn object.
data_report.add_figure(fig=montage_fig, title="GSN-HydroCel-257 Montage", caption= '256 channels')

#%% ************** Find breaks in the continuous data **********************************

break_annots = mne.preprocessing.annotate_break(RawIn, min_break_duration=20, t_start_after_previous=5, t_stop_before_next=5)
RawIn.set_annotations(RawIn.annotations + break_annots)

##*********************** Carry out low-pass and high-pass filtering ******************
hplim = .05                                                    # High-pass filter cutoff, .1Hz
Rawfilt_HP, filthp_fig, filthp_params = _dohighpass(RawIn, hplim, ch_names_all, ch_names) # Call of function to carry out high-pass filtering.

# Call of function to carry out low-pass filtering.
new_sfreq  = 250                                               # Define the sampling frequency (Hz) desired after downsampling (to be carried out after epoching)
Rawfilt_LP, filtlp_fig, filtlp_params = _dolowpass(Rawfilt_HP, new_sfreq, ch_names_all, ch_names)

data_report.add_raw (raw=RawIn, title= "Raw data", psd=False)
data_report.add_figure(fig=filthp_fig, title="High-pass filter characteristics")
data_report.add_figure(fig=filtlp_fig, title="Low-pass filter characteristics")


#%% *************** Get a summary of EOG artifacts in the data and save to the mne report. *****************************
eogs = mne.preprocessing.create_eog_epochs(Rawfilt_LP, ch_name=['E60','E19'], picks='eeg', reject_by_annotation=True).average()
eogs.apply_baseline(baseline=(None, -0.2))
fig_eog = eogs.plot_joint()
data_report.add_figure(fig=fig_eog, title='Extracted EOG Artifacts')
data_report.save(fname = report_path, overwrite=True)

#%% ***************** Carry out manual detection of noisy electrodes ************************
"""The arrow keys (up/down/left/right) can typically be used to navigate between channels and time ranges, 
    but this depends on the backend matplotlib is configured to use (e.g., mpl.use(‘TkAgg’) should work). 
    The scaling can be adjusted with - and + (or =) keys. The viewport dimensions can be adjusted with page 
    up/page down and home/end keys. Full screen mode can be to toggled with f11 key. To mark or un-mark a channel 
    as bad, click on the rather flat segments of a channel’s time series. The changes will be reflected immediately 
    in the raw object’s raw.info['bads'] entry.
"""

Events = mne.find_events(Rawfilt_LP)
mne.viz.plot_raw(Rawfilt_LP, events=Events, duration=20, n_channels=20, title='Select noisy channels', event_color='red',
                     remove_dc=True, block=True, show=True)

### ********************* Detecting Bad Electrodes *********************************************
"""
    Plot the power spectrum density (PSD) of all electrodes
"""

lensig = Rawfilt_LP.last_samp
round(len(Rawfilt_LP)/10)
pow2 = findNextPowerOf2(int(lensig/2))
chanindx = [ch_names_all.index(ic) for ic in ch_names]
print('Signal length is {} and the nearest power of two is {}'.format(lensig, pow2))
plotSpectrum(Rawfilt_LP, pow2, chanindx)
print('****plotted spectrum')

## ********************* Call of function viz_allchans() to plot all chans
wind_dur    = 10         # Duration of time interval (in seconds) presented in each window.
wind_nchans = 50         # Number of channels presented in each window.
badchans = viz_allchans(Rawfilt_LP, wind_dur, wind_nchans)
print('The channels pre-selected as bad are: {}'.format(badchans))

##%%********************* Set the average reference as projector ******************************************
"""
    Setting the average reference as a projector will permit the recalculation of the projector if we mark other 
    channels as "bad". 
"""
Rawfilt_LP.set_eeg_reference('average', projection=True)    # Create the average reference as a projector.

##%% *********************Save the filtered data with bad channels marked ********************
suffix = 'proc-filt_raw.fif'
fsave_name = '_'.join([sujname, suffix])
path2save  = os.path.join(deriv_suj_root, fsave_name)
Rawfilt_LP.save(path2save)

##%% ******************* Carry out Segmentation of the Continuous Data ************************
"""
    Note: the maximum difference between first word in sentence and critical word is 1175ms.
    So define tmin as ((1175+200)*-1)/1000 and tmax as 1.0s

"""

events = mne.find_events(Rawfilt_LP)
suj_edict   = '_'.join([sujname.split('_')[0], 'events_dict.json'])
edict_path  = os.path.join(bids_root, subdir, sesdir, suj_edict)
edict_f     = open(edict_path)
events_dict = json.load(edict_f)
max_dist = 1175/1000
tmin = (max_dist+0.2)*-1
tmax = 1.0

events_from_annots = mne.events_from_annotations(Rawfilt_LP)

events_dict2 = {int(k1): v for k1, v in events_dict.items()}
elems2delete = (1, 0, 2, 3, 4, 5, 6, 7, 8, 9)
for dki in elems2delete:
    events_dict2.pop(dki)
Epoch_data = mne.Epochs(Rawfilt_LP, events, event_id=events_dict2, tmin= tmin, tmax=tmax, baseline=(None, None), picks = chanindx, reject_by_annotation=False)


