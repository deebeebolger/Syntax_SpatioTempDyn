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
    triggers = triginfo_In["Triggers"].tolist()
    isfiller = triginfo_In["Filler?"].tolist()
    isadjadv = triginfo_In["Isadj_adv"].tolist()
    trigdata = MarksIn["trigger_code"].tolist()
    startsamp = MarksIn["start_sample"].tolist()

    P = ' '
    keyword_data = np.repeat(P, len(trigdata)).tolist()
    filler_data = np.zeros(len(trigdata), dtype=int).tolist()
    isadj_data = np.repeat(P, len(trigdata)).tolist()

    for trigcnt in range(0, len(triggers)):
        if triggers[trigcnt] in trigdata:
            tarray = np.array(trigdata)
            trigindx_curr = np.where(tarray == triggers[trigcnt])[0].tolist()
            for x in trigindx_curr:
                keyword_data[x] = keywords[trigcnt]
                filler_data[x]  = isfiller[trigcnt]
                isadj_data[x]   = isadjadv[trigcnt]
        else:
            print('current verb does not exist in current dataset')


    MarksIn['keywords'] = np.array(keyword_data)
    MarksIn['IsFiller'] = np.array(filler_data)
    MarksIn['AdjAdv']   = np.array(isadj_data)
    MarksA = MarksIn.to_numpy()

    return MarksIn, MarksA, triginfo_In

"""****************************** LOAD IN RAW EGI DATA################################"""
filename = browseFiles()
fs = filename.split('/')
datacurr = fs[-1]
basedir  = '/'.join(fs[:-1])+'/'

#**** Set the path in which to save the data and the events list*******
# These directories should be in the same folder from which you loaded the data.
Datadir  = 'Data'
Eventdir = 'EventFiles'
savedata_path   = '/'.join(fs[:-1])+'/'+Datadir+'/'
saveEvents_path = '/'.join(fs[:-1])+'/'+Eventdir+'/'

if os.path.exists(savedata_path)==FALSE:
    os.mkdir(savedata_path)
else:
    print('Directory: {} and path: {} already exists.\n' .format(Datadir, basedir ))

if os.path.exists(saveEvents_path)==FALSE:
    os.mkdir(saveEvents_path)
else:
    print('Directory: {} and path: {} already exists.\n'.format(Eventdir, basedir))


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
events_dict = dict()
Keysw       = reftrig_info["Keywords"]
Trigall     = reftrig_info["Triggers"].tolist()

for keyix in range(0,len(Keysw)):
    Kcurr = Keysw[keyix]
    kindx = np.where(Keysw == Kcurr)[0].tolist()
    Tcurr = [Trigall[i] for i in kindx]
    res   = dict({Kcurr: Tcurr})
    events_dict = {**events_dict, **res}

rawout_title = datacurr.split('.')
rawout_name  = rawout_title[0] + '.fif'
RawIn.save(savedata_path+rawout_name)

# Write dictionary to a json file
jfilepath = '/Users/bolger/PycharmProjects/Syntax_SpatioTempDyn_31-05-2022/'
jfilenom = 'events_dict.json'
json.dump(events_dict, open(jfilepath+jfilenom, 'w'))

## Write the markers_df dataframe to an excel file.
savefname1 = datacurr[:-4]  +'_EventsList.xlsx' # Change to whatever naming system you want.
with pd.ExcelWriter(saveEvents_path + savefname1) as writer:
    markers_df.to_excel(writer, sheet_name='sheet1', index=False)
