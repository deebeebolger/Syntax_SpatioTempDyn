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
savefname = 'Pilot003_EventsList.xlsx'                                     # Name of xlsx file in which to save the events list. The eprime excel file name will be generated from this.


dataIn = pd.read_excel(xls_path + fname, sheet_name='allStim')
trigsIn = pd.read_excel(xls_path + fname, sheet_name='triggers')

## Read content from dataIn dataframe
Keywords = dataIn["keyword"].tolist()
Stims = dataIn["stim"].tolist()
newID = dataIn["newID?"].tolist()
isfiller = dataIn["filler?"].tolist()
pictypes = dataIn["noun"].tolist()
trigcodes = trigsIn["triggerCode"].tolist()

Keywords_indx = list(np.arange(0, len(Keywords)))
randidx = []
wordsel = []
wordindxsel = []
IsFill = []
TrigAll = []
ImageType = []
StimType = []
FillID = []

Wsel = []
WIndxsel = []
IsFiller = []
AllTrigs = []
AllImages = []
AllStims = []
AllFillID = []


def findElements(lst1, lst2):
    return [lst1[i] for i in lst2]

tester = 0

for counter in range(0, len(Keywords)):

    print(counter)
    while tester == 0:

        np.random.seed(int(time.time()))
        randindx = np.random.choice(Keywords_indx, 1, replace=False)  # returns random indices
        currsel = findElements(Keywords, randindx)  # Find the words corresponding to current selection
        fillbool = findElements(isfiller, randindx)  # Check to see if it is a filler
        trigOrd = findElements(trigcodes, randindx)
        imageOrd = findElements(pictypes, randindx)
        stimOrd = findElements(Stims, randindx)
        fillIDord = findElements(newID, randindx)

        if counter > 0:
            isfill_pre = IsFill[counter - 1]
            if fillbool[0] == 1 and isfill_pre[0] == 1:
                print('consecutive fillers')
                tester = 0
            else:
                print('phew! No consecutive fillers')
                tester += 1

            # Need to find the last occurrence of the currsel
            if currsel in wordsel:
                x = wordsel.index(currsel)

                if (counter - x) < 5:
                    print('Oh no! The same keyword appeared in the 5 find elements of list!')
                    tester = 0
            else:
                print('Phew! No sign of this keyword in the last 5 find elements of list!')
                tester += 1

        else:
            tester += 1

    Keywords_indx.remove(randindx)

    Rindx = randindx.tolist()
    wordsel.append(currsel)
    wordindxsel.append(Rindx)
    IsFill.append(fillbool)
    TrigAll.append(trigOrd)
    ImageType.append(imageOrd)
    StimType.append(stimOrd)
    FillID.append(fillIDord)

    Wsel = np.append(Wsel, currsel)
    WIndxsel = np.append(WIndxsel, randindx)
    IsFiller = np.append(IsFiller, fillbool)
    AllTrigs = np.append(AllTrigs, trigOrd)
    AllImages = np.append(AllImages, imageOrd)
    AllStims = np.append(AllStims, stimOrd)
    AllFillID = np.append(AllFillID, fillIDord)

    tester = 0
    print(wordsel)

## ------------------------ Write the events list to a dataframe-----------------------------------------------------###
D = pd.DataFrame(
    [Wsel.tolist(), AllTrigs.tolist(), IsFiller.tolist(), AllStims.tolist(), AllImages.tolist(), AllFillID.tolist()],
    index=['Keywords', 'Triggers',
           'Filler?', 'Stims', 'ImageType', 'FillID'])
Data2Excel = D.transpose()

## ---------Search for consecutive fillers (there appears to always be one consecutive filler occurrence)------------###
searchval = [1, 1]
X = (IsFiller[:-1] == searchval[0]) & (IsFiller[1:] == searchval[1])
Xfind    = np.where(X)
allindx  = list(range(len(Data2Excel)))

## Swap with another entry if consecutive fillers have been found.
# Find the indices of the other non-repeated elements.
DExcelv2 = Data2Excel.copy()
allpos = []
newpos = []

for counter2 in range(0, len(Xfind[0])):

    np.random.seed(int(time.time()))                      # Set the random generator seed.
    poscurr = [Xfind[0][counter2], Xfind[0][counter2]+1]  # Just need to move the first one
    allpos.extend(poscurr)
    notrep   = [v for i, v in enumerate(allindx) if i not in allpos]
    tester1  = 0


    while tester1 == 0:

        alltests = []

        randfill = np.random.choice(notrep, 1, replace=False)         # Returns random indices
        print(randfill)

        # Check the random position does not have an adjacent filler
        # Check that it is not a filler position
        # Check that filler word does not appear in the 5 entries before or after this random position
        fillcurr = Data2Excel.loc[poscurr[0], 'ImageType']     # The first index of the repeated filler
        wordcurr = Data2Excel.loc[randfill[0], 'ImageType']    # The ImageType corresponding to corresponding to possible new position for poscurr
        isfill   = Data2Excel.loc[randfill[0], 'Filler?']      # Need to make sure that the new position is not a filler position.rand
        isfillpre = Data2Excel.loc[randfill[0]-1, 'Filler?']
        isfillpost = Data2Excel.loc[randfill[0]+1, 'Filler?']

        imageall = Data2Excel.loc[:, 'ImageType']
        fillindx = [idx1 for idx1 in range(len(imageall)) if imageall[idx1] == fillcurr]   # Find indices of occurrences of current filler
        wordindx = [idx2 for idx2 in range(len(imageall)) if imageall[idx2] == wordcurr]   # Find indices of occurrences of swap word.

        fdiff = np.abs(randfill - fillindx)    # Ensure that same filler word is more than 5 words apart from current filler.
        wdiff = np.abs(poscurr[0] - wordindx)  # Ensure that the same word is more than 5 words apart in the new position

        word_test  = wordcurr != fillcurr
        fill_test  = isfill == 0
        fill_pretest = isfillpre == 0
        fill_postest = isfillpost == 0
        alltests = word_test + fill_test + fill_pretest + fill_postest

        A = fdiff > 5
        A1 = wdiff > 5
        filldist1 = A.all()
        wordist1  = A1.all()

        if alltests == 4 and filldist1 and wordist1:
            print('A good candidate position')
            tester1 = 1

        elif alltests > 4 or not fill_test or not filldist1 or not wordist1:
            print('Oh no, not a good candidate position...try again!')
            tester1 = 0

    print(randfill)
    allpos.extend(randfill)

    temp  = DExcelv2.iloc[poscurr[0]].copy()
    temp2 = DExcelv2.iloc[randfill].copy()
    DExcelv2.iloc[randfill] = temp
    DExcelv2.iloc[poscurr[0]] = temp2
    randfill = []


###------Shift the fillers (Imagetype column) so that on one occurrence they correspond to auditory stimulation and on second occurrence they
## do not correspond to the auditory stimulation (Stims column)---------###

images_list = DExcelv2['ImageType'].tolist()
IsFill = DExcelv2.loc[:, 'Filler?']
IsFill = IsFill.to_list()
fillindx = [index1 for index1 in range(len(IsFill)) if IsFill[index1] == 1]  # But need to leave out the current fillers.
fill_items = [images_list[i] for i in fillindx]
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

    Fill1_imaget = Itemp1.loc[fillpair_rand[0], 'ImageType']   # Image type
    Fill1_ID     = Itemp1.loc[fillpair_rand[0], 'FillID']      # Filler ID

    Fill2_imaget = Itemp2.loc[fillrand_idx[0], 'ImageType']     # Image type
    Fill2_ID     = Itemp2.loc[fillrand_idx[0], 'FillID']        # Filler ID

    DExcelv2.loc[fillpair_rand[0], 'ImageType'] = Fill2_imaget
    DExcelv2.loc[fillpair_rand[0], 'FillID']    = Fill2_ID

    DExcelv2.loc[fillrand_idx, 'ImageType'] = Fill1_imaget
    DExcelv2.loc[fillrand_idx, 'FillID']    = Fill1_ID

    # The current filler word and the swapped filler word are excluding from filler word lists.
    allfills.append(fillitems_u[fcount])
    allfills.append(fword)

    # Redefine the fill item unique list
    fillitems_u = [fillcnt for fillcnt in fillitems_u if fillcnt not in allfills]

    if len(fillitems_u) == 0:
        tester = 1
    else:
        tester = 0


###----------------------Write the Dataframe to an Excel file-------------------------------------------------------####

#with pd.ExcelWriter(xls_path + savefname) as writer:
#    DExcelv2.to_excel(writer, sheet_name='sheet1')

### Resume the inter-keyword interval by extracting the indices of occurrence of each keyword.###
AllKeywords_intval = []
for wcnt in range(0, len(wordsel)):
    Windices = [index for (index, item) in enumerate(wordsel) if item == wordsel[wcnt]]
    wdiff = Windices[1] - Windices[0]
    print("Inter-keyword interval for current word %s is %d" % (wordsel[wcnt], wdiff))
    AllKeywords_intval.append(wdiff)


###-------------------Restructure the Dataframe to Eprime format and save to Excel file-----------------------------####
# The Excel file should have the following headers: ID; Weight (1); Nested; Procedure (trialProc); Audio; picture;
# triggerAudio; triggerFixation

# Prepare the audio column. It needs to show path and .*wav extension.
allfillers    = DExcelv2['Filler?'].tolist()
fillindex1    = [fillindx1 for fillindx1 in range(len(allfillers)) if allfillers[fillindx1] == 1]   # Find the indiices of the fillers
audioID = DExcelv2['FillID']
audio   = 'a/'+audioID+'.wav'
for acnt in range(0, len(fillindex1)):
    audio[fillindex1[acnt]] = 'a/'+audioID[fillindex1[acnt]]+'_F.wav'

# Prepare the ID column.
IDcol = np.linspace(1,len(audio), len(audio))

# Prepare the Weight column.
Weightcol = np.ones(len(audio), dtype=int)

# Prepare the Procedure column.
P = 'trialProc'
Proccol = np.repeat(P, len(audio))

# Prepare the Picture column
imtype = DExcelv2['ImageType'].tolist()
qmark  = '?'
piccol = np.repeat(qmark, len(audio)).T.tolist()
for pcnt in range(0, len(fillindex1)):
    piccol[fillindex1[pcnt]] = 'pictures/'+imtype[fillindex1[pcnt]]+'.png'

# To prepare the "correct" column test if "ImageType" is in "Stims":
# If yes ==> correct (Y), if no ==> incorrect (N).
stimtype = DExcelv2['Stims'].tolist()
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
Deprime.loc[48.5] = [49, 1, '?', 'breakProc', '?', '?', '?', '?', '?']
Deprime.loc[96.5] = [97, 1, '?', 'breakProc', '?', '?', '?', '?', '?']
Deprime = Deprime.sort_index().reset_index(drop=True)
IDcol2 = np.linspace(1,len(audio)+2, len(audio)+2)
IDcol2 = IDcol2.T
Deprime['ID'] = IDcol2

savefname1 = savefname[:-5]  +'_eprime.xlsx' # Change to whatever naming system you want.
with pd.ExcelWriter(xls_path + savefname1) as writer:
    Deprime.to_excel(writer, sheet_name='sheet1', index=False)
