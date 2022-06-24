# SCRIPT TO GENERATE EVENT LISTS FOR EXPERIMENTAL PROTOCOL.
"""
    Programmed by: Deirdre BOLGER
    Date: May 2022
    The events lists must adhere to the following rules:
    - No adjacent fillers
    - No stimuli with the same keyword less than 5 trials apart.
    For each participant, there are two blocks (144*2 trials)
    - The fillers are assigned such that on one presentation they correspond to the audio
    stimulation and on a second presentation they do not.
    Requires: the xlsx file "Experimental_Lists_Bissera1.xlsx"
    Output:
    The script outputs two excel files:
    - An excel file (xlsx) with events ordered according to above rules. But this has been commented out.
    - An excel file (xlsx) with Eprime structure to be imported into Eprime script.

    Important: sometimes the script will not manage to converge on the correct organisation of events and will continue to search
    and search. In this case, just stop the script running and try to run it again. Around 99% of the time it works.

"""

import pandas as pd
import numpy as np
import time
import openpyxl
from array import array

# Define filenames and paths
fname = 'Experimental_Lists_Bissera1.xlsx'
xls_path = r'/Users/bolger/Documents/work/Projects/SpatioTempDyn_Syntax/'  # Path in which to save output xlsx files and from which to load fname.
savefname = 'Pilot0099_EventsList.xlsx'                                     # Name of xlsx file in which to save the events list. The eprime excel file name will be generated from this.

#dataIn = pd.read_excel(xls_path + fname, sheet_name='allStim')
trigsIn = pd.read_excel(xls_path + fname, sheet_name='triggers')

## Read content from dataIn dataframe
Keywords   = trigsIn["keyword"].tolist()
StimID     = trigsIn["stimID"].tolist()
Conds      = trigsIn["condition"].tolist() #Fillers = 3 (adj, adv), adj = 1, adv = 2
nouns      = trigsIn["noun"].tolist()
stims_sent = trigsIn["stim"].tolist()   # the auditory phrase presented
trigcodes  = trigsIn["triggerCode"].tolist()

StimID_indx = list(np.arange(0, len(StimID)))

Rindx = []         # Randomly chosen indices
StimIndex_sel = [] # Randomly chosen entry of
KeyWord_sel = []
Cond_sel = []
AudStim_sel = []
Trig_sel = []
Noun_sel = []
Phrase_sel = []

wordindxsel = []
IsFill = []
TrigAll = []
ImageType = []
StimType = []
FillID = []
FillID2 = []

Wsel = []
WIndxsel = []
IsFiller = []
AllTrigs = []
AllImages = []
AllStims = []
AllFillID = []
AllID = []

#def findElements(lst1, lst2):
#    return [lst1[i] for i in lst2]

def test_keywords(Allkeywords, CurrKeywords, cntr, test1):
    if CurrKeywords in Allkeywords:
        print('Looking for neighbouring keywords')
        x = Allkeywords.index(keywordsel)

        if (cntr - x) < 5:
            print('Oh no! The same keyword appeared in the 5 find elements of list!')
            test2 = 0
        else:
            print('Phew! No sign of this keyword in the last 5 find elements of list!')
            test2 = 1

        test = test1 * test2
    else:
        test = test1
    return test



tester = 0
for counter in range(0, len(StimID)):

    print(counter)
    while tester == 0:

        currsel = []
        p = np.ones(len(StimID_indx), dtype=int)
        p = p/len(StimID_indx)
        randindx = np.random.choice(StimID_indx, 1, replace=False, p = p)  # Returns random indices
        currsel  = StimID_indx.index(randindx[0])             # Find the word index corresponding to current random selection
        audstim  = StimID[currsel]                       # Find the currently selected audio stimulus.
        condsel  = Conds[currsel]                   # Find the currently selected filler type (0 or 1)
        trigsel  = trigcodes[currsel]                    # Find the currently selected trigger code.
        stimsel  = stims_sent[currsel]               # Find the currently selected audio sentence.
        keywordsel = Keywords[currsel]               # Find the currently selected keyword
        nounsel    = nouns[currsel]                    # Find the currently selected noun (picture to be shown on filler trials)

        if counter > 0:
            isfill_pre = Cond_sel[counter - 1]

            if condsel == 3 and isfill_pre == 3:
                print('consecutive fillers')
                tester1 = 0

                tester = test_keywords(KeyWord_sel, keywordsel, counter, tester1)

            else:
                print('phew! No consecutive fillers')
                tester1 = 1
                tester = test_keywords(KeyWord_sel, keywordsel, counter, tester1)

        else:

            tester += 1

    StimID_indx.remove(randindx[0])  # Remove the currently selected entry index.
    StimID.pop(currsel)           # Remove the currently selected and assigned audio stimulus
    Conds.pop(currsel)   # Remove the currently selected and assigned filler
    trigcodes.pop(currsel)        # Remove the currently selected and assigned trigger code.
    stims_sent.pop(currsel)       # Remove the currently selected and assigned sentence stimulus
    Keywords.pop(currsel)      # Remove the currently selected and assigned keyword
    nouns.pop(currsel)            # Remove the currently selected and assigned noun.

    Rindx.append(randindx[0])             # Add the randomly selected index list
    StimIndex_sel.append(randindx[0])  # Add to the randomly assigned stimulus index list.
    KeyWord_sel.append(keywordsel)     # Add the the randomly assigned keyword list
    Cond_sel.append(condsel)         # Add to the randomly assigned filler list
    AudStim_sel.append(audstim)        # Add to the randomly assigned auditory stimulus list
    Trig_sel.append(trigsel)       # Add to the randomly assigned trigger code list
    Phrase_sel.append(stimsel)  # Add to the randomly assigned audio stimulus phrase list.
    Noun_sel.append(nounsel)       # Add to the randomly assigned nouns list.

    tester = 0
    print(AudStim_sel)

## ------------------------ Write the events list to a dataframe-----------------------------------------------------###
D = pd.DataFrame(
    [KeyWord_sel, Trig_sel, Cond_sel, AudStim_sel, Noun_sel, Phrase_sel, Rindx],
    index=['Keywords', 'Triggers',
           'Filler?', 'StimAudio', 'Picture', 'StimulusPhrase', 'RandIndx'])
Data2Excel = D.transpose()

## ---------Search for consecutive fillers (there appears to always be one consecutive filler occurrence)------------###
searchval = [3, 3]
X = (IsFiller[:-1] == searchval[0]) & (IsFiller[1:] == searchval[1])
Xfind    = np.where(X)


###------Shift the fillers (Pictures column) so that on one occurrence they correspond to auditory stimulation and on second occurrence they
## do not correspond to the auditory stimulation (Stims column)---------###
## ONLY NEED TO CHANGE THE PICTURE COLUMN

DExcelv2 = Data2Excel
pics_list = DExcelv2['Picture'].tolist()
IsFill = DExcelv2.loc[:, 'Filler?']
IsFill = IsFill.to_list()
fillindx = [index1 for index1 in range(len(IsFill)) if IsFill[index1] == 3]  # But need to leave out the current fillers.
fill_items = [pics_list[i] for i in fillindx]

# Correct the sorciere word.
iwitch = [iw for iw in range(0, len(fill_items)) if fill_items[iw] =='sorcière']
for iw in iwitch:
    fill_items[iw] = 'sorciere'

fillitems_u = np.unique(fill_items)
fillitems_u = fillitems_u.tolist()
allfills = []

fcount = 0
tester = 0

while tester == 0:

    print(fillitems_u[fcount])      # Print to screen the current

    indexes = [index for index in range(
        len(fill_items)) if fill_items[index] == fillitems_u[fcount]]    # Find the indices of the  2 occurrences of each filler picture.
    I = [fillindx[i1] for i1 in indexes]
    print(I)

    # Randomly select one of the filler pair (I) to become non-compatible.
    np.random.seed(int(time.time()))
    fillpair_rand = np.random.choice(I, 1, replace=False)

    # Randomly select one of the other fillers to exchange with fillpair_rand
    # Find all the other pairs other than the current filler pair (I)
    fillother_indx = [e for e in fillindx if e not in I]
    fillrand_idx = np.random.choice(fillother_indx, 1, replace=False)   # The randomly chosen index of other filler for exchange
    findx = fillindx.index(fillrand_idx)                                # The filler index of the randomly chosen filler
    fword = fill_items[findx]                                           # The randomly chosen filler word

    # swap between positions fillpair_rand and fillrand_idx
    Itemp1 = DExcelv2.iloc[fillpair_rand]
    Itemp2 = DExcelv2.iloc[fillrand_idx]

    Fill1_imaget = Itemp1.loc[fillpair_rand[0], 'Picture']   # Image type
    Fill2_imaget = Itemp2.loc[fillrand_idx[0], 'Picture']     # Image type

    if Fill1_imaget == 'sorcière':
        Fill1_imaget = 'sorciere'

    if Fill2_imaget == 'sorcière':
        Fill2_imaget = 'sorciere'

    DExcelv2.loc[fillpair_rand[0], 'Picture'] = Fill2_imaget
    DExcelv2.loc[fillrand_idx, 'Picture'] = Fill1_imaget

    # The current filler word and the swapped filler word are excluded from filler word lists.
    allfills.append(fillitems_u[fcount])
    allfills.append(fword)

    # Redefine the fill item unique list
    fillitems_u = [fillcnt for fillcnt in fillitems_u if fillcnt not in allfills]

    if len(fillitems_u) == 0:
        tester = 1
    else:
        tester = 0


###----------------------Write the Dataframe to an Excel file-------------------------------------------------------####

with pd.ExcelWriter(xls_path + savefname) as writer:
    DExcelv2.to_excel(writer, sheet_name='sheet1')

### Resume the inter-keyword interval by extracting the indices of occurrence of each keyword.###
AllKeywords_intval = []
for wcnt in range(0, len(KeyWord_sel)):
    Windices = [index for (index, item) in enumerate(KeyWord_sel) if item == KeyWord_sel[wcnt]]
    wdiff = Windices[1] - Windices[0]
    print("Inter-keyword interval for current word %s is %d" % (KeyWord_sel[wcnt], wdiff))
    AllKeywords_intval.append(wdiff)


###-------------------Restructure the Dataframe to Eprime format and save to Excel file-----------------------------####
# The Excel file should have the following headers: ID; Weight (1); Nested; Procedure (trialProc); Audio; picture;
# triggerAudio; triggerFixation

# Prepare the audio column. It needs to show path and .*wav extension.
allfillers    = DExcelv2['Filler?'].tolist()
fillindex1    = [fillindx1 for fillindx1 in range(len(allfillers)) if allfillers[fillindx1] == 3]   # Find the indiices of the fillers
audioID = DExcelv2['StimAudio']
audio   = 'a/'+audioID+'.wav'

# Prepare the ID column.
IDcol = np.linspace(1,len(audio), len(audio))

# Prepare the Weight column.
Weightcol = np.ones(len(audio), dtype=int)

# Prepare the Procedure column.
P = 'fillerProc'
Proccol = np.repeat(P, len(audio))
not_fillindex1 = [fillindx1 for fillindx1 in range(len(allfillers)) if allfillers[fillindx1] < 3]
for pccnt in range(0, len(not_fillindex1)):
    Proccol[not_fillindex1[pccnt]] = 'trialProc'

# Prepare the Picture column
imtype = DExcelv2['Picture'].tolist()
qmark  = '?'
piccol = np.repeat(qmark, len(audio)).T.tolist()
for pcnt in range(0, len(fillindex1)):
    piccol[fillindex1[pcnt]] = 'pictures/'+imtype[fillindex1[pcnt]]+'.png'

# To prepare the "correct" column test if "ImageType" is in "Stims":
# If yes ==> correct (Y), if no ==> incorrect (N).
stimtype = DExcelv2['StimulusPhrase'].tolist()
corrcol = np.repeat(qmark, len(audio)).T.tolist()
for corrcnt in range(0, len(fillindex1)):
    if imtype[fillindex1[corrcnt]] in stimtype[fillindex1[corrcnt]]:
        print('Correct!')
        corrcol[fillindex1[corrcnt]] = 'Y'
    else:
        print('Incorrect')
        corrcol[fillindex1[corrcnt]] = 'N'

TriggerCol = DExcelv2['Triggers']
trigfixCol = np.ones(len(audio), dtype=int)
nestedCol = np.repeat('', len(audio)).tolist()

# Create a dataframe to write to Eprime.
Deprime = pd.DataFrame(
    [IDcol.T.tolist(), Weightcol.T.tolist(), nestedCol, Proccol.T.tolist(), audio.T.tolist(), piccol, corrcol, TriggerCol.T.tolist(), trigfixCol.T.tolist()],
    index=['ID', 'Weight', 'Nested', 'Procedure', 'Audio', 'picture', 'correct', 'triggerAudio', 'triggerFixation'])
Deprime = Deprime.T   # Transpose of the dataframe

# Add in a  trialProc to breakProc on sample 48+1 and 96+1...
Deprime.loc[48.5] = [49, 1, '?', 'breakProc', '?', '?', '?', '200', '?']
Deprime.loc[96.5] = [97, 1, '?', 'breakProc', '?', '?', '?', '200', '?']
Deprime = Deprime.sort_index().reset_index(drop=True)
IDcol2 = np.linspace(1,len(audio)+2, len(audio)+2)
IDcol2 = IDcol2.T
Deprime['ID'] = IDcol2

savefname1 = savefname[:-5]  +'_eprime.xlsx' # Change to whatever naming system you want.
with pd.ExcelWriter(xls_path + savefname1) as writer:
    Deprime.to_excel(writer, sheet_name='sheet1', index=False)
