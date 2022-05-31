# SCRIPT TO CORRECT THE TRIGGERS
"""
    Programmed by: Deirdre BOLGER
    Date: 31 May 2022
    Load in the excel file with the triggers (both old and the new).

"""

import pandas as pd
import numpy as np
import openpyxl
import xlwt


fname = 'experimental_lists_Bissera.xlsx'
xls_path = r'/Users/bolger/Documents/work/Projects/SpatioTempDyn_Syntax/'

trigsIn    = pd.read_excel(xls_path + fname, sheet_name='triggers')        # Load in the triggers and StimID variables.
newtrigsIn = pd.read_excel(xls_path + fname, sheet_name='new_triggers')    # Load in the triggerAudio and Audio columns.

# EXTRACT THE DATA COLUMNS FROM EXCEL SHEETS (TRIGGERS) AND (NEW_TRIGGERS)
old_trigs = trigsIn['triggerCode'].tolist()
oldtrig_ID = trigsIn['stimID'].tolist()
new_trigs = newtrigsIn['triggerAudio'].tolist()
newtrig_ID = newtrigsIn['Audio'].tolist()

# Compare the oldtrig_ID and newtrig_ID columns.
oldtrigs_bis = old_trigs
for counter in range(0, len(newtrig_ID)):

    trigtesta = newtrig_ID[counter]
    trigtest = trigtesta[2:-4]
    if trigtest[-1] == 'F':
        trigtest = trigtesta[2:-6]
    match_indx = oldtrig_ID.index(trigtest)

    oldtrigs_bis[match_indx] = new_trigs[counter]

## APPEND THE NEW TRIGGER CODES TO THE CURRENT EXCEL FILE (TRIGGERS SHEET)

TriggersNew = pd.DataFrame(oldtrigs_bis)
writer = pd.ExcelWriter(xls_path + fname)

TriggersNew.to_excel(writer, sheet_name='NewTrigs')
writer.save()

