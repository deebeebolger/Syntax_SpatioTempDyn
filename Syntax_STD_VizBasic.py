from tkinter import filedialog
import os
import mne
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
mne.set_config('MNE_BROWSER_BACKEND','matplotlib')


def browseFiles():
    cdpath = os.getcwd()
    filename = filedialog.askopenfilename(initialdir=cdpath, title="Select a  file",
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
    mne.viz.plot_raw(raw_eeg, duration=vizdur, n_channels=channum_viz, title='Raw EEG Data', event_color='red',
                     remove_dc=True, block=True, show=True)

    bads = raw_eeg.info['bads']
    return bads


###*************************LOAD IN THE FILE MANUALLY
filename = browseFiles()
dataIn = mne.read_epochs(filename, preload=True, verbose=None)
event_id = dataIn.event_id
dataIn.plot(block=True)


###************************Load in continuous data manually***********************
filename_cont = browseFiles()
dataIn_cont = mne.io.read_raw_fif(filename_cont, preload=True, verbose=None)
wind_dur    = 10         # Duration of time interval (in seconds) presented in each window.
wind_nchans = 50         # Number of channels presented in each window.
badchans = viz_allchans(dataIn_cont, wind_dur, wind_nchans)
print('The channels pre-selected as bad are: {}'.format(badchans))

## ********************* Call of function viz_allchans() to plot all chans****************
wind_dur    = 60         # Duration of time interval (in seconds) presented in each window.
wind_nchans = 50         # Number of channels presented in each window.
badchans = viz_allchans(dataIn, wind_dur, wind_nchans)
