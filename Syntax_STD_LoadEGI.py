import mne
import numpy as np
from datetime import datetime
from glob import glob
from os.path import basename, join, splitext
from xml.etree.ElementTree import parse
import pandas as pd
import json as json
from tkinter import *
from tkinter import filedialog
import os

def browseFiles():
    cdpath = os.getcwd()
    filename = filedialog.askopenfilename(initialdir = cdpath, title="Select a raw file", filetypes=[("all files", "*.*")])
    return filename

def _parse_xml(xml_file):
    """Parse XML file."""
    xml = parse(xml_file)
    root = xml.getroot()
    return _xml2list(root)

def _xml2list(root):
    """Parse XML item."""
    output = []
    for element in root:

        if len(element) > 0:
            if element[0].tag != element[-1].tag:
                output.append(_xml2dict(element))
            else:
                output.append(_xml2list(element))

        elif element.text:
            text = element.text.strip()
            if text:
                tag = _ns(element.tag)
                output.append({tag: text})
    return output

def _ns(s):
    """Remove namespace, but only if there is a namespace to begin with."""
    if '}' in s:
        return '}'.join(s.split('}')[1:])
    else:
        return s


def _xml2dict(root):
    """Use functions instead of Class.

    remove namespace based on
    http://stackoverflow.com/questions/2148119
    """
    output = {}
    if root.items():
        output.update(dict(root.items()))

    for element in root:
        if len(element) > 0:
            if len(element) == 1 or element[0].tag != element[1].tag:
                one_dict = _xml2dict(element)
            else:
                one_dict = {_ns(element[0].tag): _xml2list(element)}

            if element.items():
                one_dict.update(dict(element.items()))
            output.update({_ns(element.tag): one_dict})

        elif element.items():
            output.update({_ns(element.tag): dict(element.items())})

        else:
            output.update({_ns(element.tag): element.text})
    return output

def _ns2py_time(nstime):
    """Parse times."""
    nsdate = nstime[0:10]
    nstime0 = nstime[11:26]
    nstime00 = nsdate + " " + nstime0
    pytime = datetime.strptime(nstime00, '%Y-%m-%d %H:%M:%S.%f')
    return pytime

def _combine_triggers(data, remapping=None):
    """Combine binary triggers."""
    new_trigger = np.zeros(data.shape[1])
    if data.astype(bool).sum(axis=0).max() > 1:  # ensure no overlaps
        logger.info('    Found multiple events at the same time '
                    'sample. Cannot create trigger channel.')
        return
    if remapping is None:
        remapping = np.arange(data) + 1
    for d, event_id in zip(data, remapping):
        idx = d.nonzero()
        if np.any(idx):
            new_trigger[idx] += event_id
    return new_trigger

def load_triginfo(pathin, fnamein, MarksIn):
    """Load in the excel resuming the trigger information"""
    triginfo_In = pd.read_excel(pathin + fnamein, sheet_name=0)

    keywords = triginfo_In["Keywords"].tolist()
    stimID   = triginfo_In["StimID"].tolist()
    triggers = triginfo_In["Triggers"].tolist()
    isfiller = triginfo_In["Filler?"].tolist()
    isadjadv = triginfo_In["Isadj_adv"].tolist()
    trigdata = MarksIn["trigger_code"].tolist()
    startsamp = MarksIn["start_sample"].tolist()

    P = ' '
    keyword_data = np.repeat(P, len(trigdata)).tolist()
    stimID_data  = np.repeat(P, len(trigdata)).tolist()
    filler_data = np.zeros(len(trigdata), dtype=int).tolist()
    isadj_data = np.repeat(P, len(trigdata)).tolist()

    for trigcnt in range(0, len(triggers)):
        if triggers[trigcnt] in trigdata:
            tarray = np.array(trigdata)
            trigindx_curr = np.where(tarray == triggers[trigcnt])[0].tolist()
            for x in trigindx_curr:
                keyword_data[x] = keywords[trigcnt]
                stimID_data[x]  = stimID[trigcnt]
                filler_data[x]  = isfiller[trigcnt]
                isadj_data[x]   = isadjadv[trigcnt]
        else:
            print('current verb does not exist in current dataset')


    MarksIn['keywords'] = np.array(keyword_data)
    MarksIn['stimID']   = np.array(stimID_data)
    MarksIn['IsFiller'] = np.array(filler_data)
    MarksIn['AdjAdv']   = np.array(isadj_data)
    MarksA = MarksIn.to_numpy()

    return MarksIn, MarksA, triginfo_In

def create_bidsevents(markersIn, savepath, sujname_curr):

    # Create a dataframe from the following lists.
    onsets = markersIn['start'].T.tolist()
    durcol = np.repeat('n/a', len(onsets)).tolist()
    samps  = markersIn['start_sample'].T.tolist()
    trigs  = markersIn['trigger_code'].T.tolist()
    keywd  = markersIn['keywords'].T.tolist()
    fillr  = markersIn['IsFiller'].T.tolist()
    wtype  = markersIn['AdjAdv'].T.tolist()

    allists = [onsets, durcol, samps, keywd, fillr, wtype, trigs]
    cols    = ["onsets", "duration", "sample", "keyword", "IsFiller", "Word_type", "trigger_codes"]
    eventsDF = pd.DataFrame(np.transpose(allists), columns = ["onsets", "duration", "sample", "keyword", "IsFiller", "Word_type", "trigger_codes"])
    print(eventsDF)
    currnom = 'sub-' + sujname_curr[:-4] + '_events.tsv'
    savepath_curr = savepath + currnom
    eventsDF.to_csv(savepath_curr, sep="\t")

    return eventsDF, currnom


"""****************************** LOAD IN RAW EGI DATA################################"""
filename = browseFiles()
fs = filename.split('/')
datacurr = fs[-1]
basedir  = '/'.join(fs[:-2])
print('The base directory is: {}\n' .format(basedir))

# Load in the onsets_file
fonsets = 'onset_times.xlsx'
onset_path = '/'.join([basedir, fonsets])
Onset_info = pd.read_excel(onset_path)
Onset_stim = Onset_info['STIM']
Onset_adjst = Onset_info['ADJ_START']
Onset_advst = Onset_info['ADV_START']

#**** Set the path in which to save the data and the events list*******
# These directories should be in the same folder from which you loaded the data.
Datadir  = 'Data'
Eventdir = 'EventFiles'
savedata_path   = basedir +'/' + Datadir +'/'
saveEvents_path = basedir +'/' + Eventdir +'/'

if os.path.exists(savedata_path)==FALSE:
    os.mkdir(savedata_path)
else:
    print('Directory: {} and path: {} already exists.\n' .format(Datadir, basedir ))

if os.path.exists(saveEvents_path)==FALSE:
    os.mkdir(saveEvents_path)
else:
    print('Directory: {} and path: {} already exists.\n'.format(Eventdir, basedir))

##****** CREATE THE STUDY DIRECTORY*****
study_name = 'SyntaxDyn_hd'
sdir       = [basedir, study_name]
studydir   =  '/'.join(sdir)

# Create a directory for the current participant.
if datacurr[1] =='_':
    newsub_dir = 'sub-' + '0' + datacurr[0]
else:
    newsub_dir = 'sub-' + datacurr[0:2]
save_currsuj_path = '/'.join([studydir, newsub_dir])
if os.path.exists(save_currsuj_path)==FALSE:
    os.mkdir(save_currsuj_path)
else:
    print('Directory: {} and path: {} already exists.\n'.format(newsub_dir, save_currsuj_path))

## Create a session directory within the newsub_dir directory
session_dir = 'ses-1'
save_sessionpath = '/'.join([save_currsuj_path, session_dir])
if os.path.exists(save_sessionpath)==FALSE:
    os.mkdir(save_sessionpath)
else:
    print('Directory: {} and path: {} already exists.\n'.format(session_dir, save_currsuj_path))

# Create an "eeg" folder within the current directory in which to save eeg data.
eeg_dir = 'eeg'
saveeeg_dir = '/'.join([save_sessionpath, eeg_dir])
if os.path.exists(saveeeg_dir)==FALSE:
    os.mkdir(saveeeg_dir)
else:
    print('Directory: {} and path: {} already exists.\n'.format('eeg', saveeeg_dir))

#*********Set the search path for the TriggerCoding_Summary.xlsx file to the current script directory******##
trigxl_path_all = os.path.abspath('Syntax_STD_LoadEGI.py')
trigXL          = trigxl_path_all.split('/')
trigxl_path     = '/'.join(trigXL[:-1])+'/'
trigxl_file     = 'TriggerCoding_SummaryJuly2022.xlsx'

RawIn = mne.io.read_raw_egi(filename, channel_naming='E%d', verbose=None, preload=True)   # Load in raw EGI data in *.mff format
sfreq = RawIn.info['sfreq']   # get the sampling frequency

"""Extract the events.

Parameters
----------
filename : str
    File path.
sfreq : float
    The sampling frequency
"""
orig = {}
for xml_file in glob(join(filename, '*.xml')):     # Extracting the xml files composing the current *.mff file.
    xml_type = splitext(basename(xml_file))[0]
    orig[xml_type] = _parse_xml(xml_file)
xml_files = orig.keys()
xml_events = [x for x in xml_files if x[:7] == 'Events_']
for item in orig['info']:
    if 'recordTime' in item:
        start_time = _ns2py_time(item['recordTime'])
        break
markers = []
code = []
for xml in xml_events:
    for event in orig[xml][2:]:
        event_start = _ns2py_time(event['beginTime'])
        start = (event_start - start_time).total_seconds()
        if event['code'] not in code:
            cc = event['code'][1:]
            x = sum(ccint.isalpha() for ccint in cc)
            if x == 0:
                code_temp = int(event['code'][1:])
                curr_code = 255 - code_temp
                code.append(curr_code)
            elif x == 1:
                code_temp = int(event['code'][2:])
                curr_code = 255 - code_temp
                code.append(curr_code)
            elif x > 1:
                curr_code = 999
                code.append(curr_code)
        marker = {'DINs': event['code'],
                  'trigger_code': curr_code,
                  'start': start,
                  'start_sample': int(np.fix(start * sfreq)),
                  'end': start + float(event['duration']) / 1e9,
                  'chan': None,
                  }
        markers.append(marker)  # Contains information regarding all event markers (in dict format).
events_tims = dict()
for ev in code:
    trig_samp = list(c['start_sample'] for n,
                     c in enumerate(markers) if c['trigger_code'] == ev)
    events_tims.update({ev: trig_samp})    # Defines the onset times of each event in dict format.

markers_df = pd.DataFrame.from_dict(markers)

# Call of function to add trigger information to the dataframe.
markers_df, markers_nda, reftrig_info = load_triginfo(trigxl_path, trigxl_file, markers_df)

# Add the corrected markers to the raw object.
E = mne.find_events(RawIn)
E[:,2] = markers_nda[:,1]
RawIn.add_events(E, 'STI 014', replace=True)

# Create an events dictionary
events_dict  = dict()
Keysw_notlist = reftrig_info["Keywords"]
Keysw        = reftrig_info["Keywords"]
Trigall      = reftrig_info["Triggers"].tolist()
IsAdjAdv     = reftrig_info["Isadj_adv"].tolist()
IsFiller_all = reftrig_info["Filler?"].tolist()
AdjAdv_wd  = reftrig_info["Adj_Advb"].tolist()
StimID_all = reftrig_info["StimID"]

for keyix in range(0,len(Keysw)):
    Kcurr        = Keysw[keyix]
    kindx        = np.where( Keysw_notlist == Kcurr)[0].tolist()
    StID_curr    = [StimID_all[i1] for i1 in kindx]
    StimID_curr  = StID_curr[0].split('_')[0]

    oindx = np.where(Onset_stim == StimID_curr)[0].tolist()  # adding onset of adverbs and adjectives to dict
    Onset_adj = Onset_adjst[oindx]
    Onset_adv = Onset_advst[oindx]

    Tcurr        = [Trigall[i] for i in kindx]
    AdjAdv_curr  = [IsAdjAdv[i2] for i2 in kindx]
    AdjAdv_words = [AdjAdv_wd[i3] for i3 in kindx]
    Filler_curr  = [IsFiller_all[i4] for i4 in kindx]
    Keyword_curr  = [Keysw[i5] for i5 in kindx]
    res   = dict({StimID_curr: [Keyword_curr, Tcurr, AdjAdv_curr, AdjAdv_words, Filler_curr,
                                [Onset_adj.values[0].tolist(), Onset_adv.values[0].tolist()]]})
    events_dict = {**events_dict, **res}   ## Could be used as a mapping

## Create a mapping for the events to annotations procedure
events_mapping = dict()
Trig_askey = reftrig_info["Triggers"]
Kwords_all = reftrig_info["Keywords"].tolist()

for trigix in range(0, len(Trig_askey)):
    Trgcurr  = Trig_askey[trigix]
    res1     = dict({Trgcurr: Kwords_all[trigix]})
    events_mapping = {**events_mapping, **res1}

# Need a call of function to create the events.tsv file for the current participant.
# note that this file should be saved with the dataset (.fif)
eventbids, evt_title = create_bidsevents(markers_df, save_sessionpath, datacurr)
eventsall = eventbids.values
annots_from_events = mne.annotations_from_events(events = E, event_desc=events_mapping, sfreq=RawIn.info['sfreq'], orig_time=RawIn.info['meas_date'])

# Store any annotations already present concerning empty data intervals ('BAD_ACQ_SKIP")
annot_orig = RawIn.annotations
RawIn.set_annotations(annots_from_events + annot_orig)

## Save the current subject raw file
rawout_title = datacurr.split('.')
rawout_title2 = rawout_title[0].split('_')
if len(rawout_title2[0]) ==1:
    rawout_corr = '0'+rawout_title2[0]
else:
    rawout_corr = rawout_title2[0]
rawout_name  = 'sub-'+rawout_corr + '_' + session_dir + '_' + eeg_dir + '.fif'
RawIn.save('/'.join([saveeeg_dir,rawout_name]), overwrite=TRUE)

# Write dictionary to a json file
jfilepath = save_sessionpath + '/'
jfilenom = 'sub-'+ rawout_corr + '_events_dict.json'
json.dump(events_dict, open(jfilepath+jfilenom, 'w'))

## Write the markers_df dataframe to an excel file.
savefname1 = datacurr[:-4]  +'_EventsList.xlsx' # Change to whatever naming system you want.
with pd.ExcelWriter(saveEvents_path + savefname1) as writer:
    markers_df.to_excel(writer, sheet_name='sheet1', index=False)
