import mne

fpath = '/Users/bolger/Documents/work/Projects/SpatioTempDyn_Syntax/Data/'
datacurr = '013_20220608_030949.mff'
fullpath = fpath + datacurr

RawIn = mne.io.read_raw_egi(fullpath, channel_naming='E%d', verbose=None)   # Do not need to specify the
events_all = mne.find_events(RawIn, stim_channel='STI 014', output='onset', consecutive=False)
mne.viz.plot_events(events_all)

RawIn.plot(events_all, color='gray')