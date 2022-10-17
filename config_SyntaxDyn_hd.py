"""
"""
import os
import mne
import mne_bids
from mne_bids import BIDSPath, read_raw_bids
import json
import numpy as np

study_name = 'DynSyn_EEG_PC'
study_root = '/Users/bolger/Documents/work/Projects/Project_DynSyn_EEG_PC'
bids_root  =  '/'.join([study_root, study_name])    # Root directory
deriv_root = '/'.join([bids_root, 'derivatives'])

sessions = ['1']
subjects = ['10']
ch_types = ['eeg']
data_type = 'eeg'

interactive = False
resample_sfreq = 250
eeg_template_montage = mne.channels.make_standard_montage('GSN-HydroCel-257')

l_freq = 0.1
h_freq = 40

find_breaks = False
min_break_duration = 4
t_break_annot_start_after_previous_event = 2.0
t_break_annot_stop_before_next_event = 1.5

ica_reject = None #dict(eeg=350e-6, eog=500e-6)
reject = 'autoreject_global'
eeg_reference = 'average'

random_state = 42
spatial_filter = None
ica_max_iterations = 1000
ica_eog_threshold = 2

run_source_estimation = False

on_error = 'abort'
on_rename_missing_events = 'warn'

analyze_channels = 'all'
plot_psd_for_runs = []

N_JOBS = 2
parallel_backend = 'loky'

print('got this far')

## READ IN THE DICTIONAY OF EVENT INFORMATION STORED AS A JSON FILE.
#  We load in the events_dict.son file of the first subject in the list as the
#  events_dict.json is the same for each subject.
#  The events_dict.json file is stored at the level of 'ses-1' directory.

trigcorr_fpath = '/'.join([bids_root, 'sub-'+subjects[0], 'ses-1'])
onsetfile = [f for f in os.listdir(trigcorr_fpath) if f.endswith('dict.json')]
with open('/'.join([trigcorr_fpath, onsetfile[0]])) as fonsets:
    onsetcorr_data = fonsets.read()
onsetcorr_js = json.loads(onsetcorr_data)   # Reconstruct onsetcorr_data as a dictionary.

K = list(onsetcorr_js.keys())
allconds  = [onsetcorr_js[key] for key in K]
epochs_tim = -0.25
epochs_tmax = 1.000
baseline = (None, 0)
conditions = allconds
contrasts = [('rouge/adj', 'rouge/adv')]
# report_evoked_n_time_points = 20
