#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import os
from datetime import datetime

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('pridedb_version', metavar='pridedb_version', type=str, help="Path to raw folder with the data.")


# In[]:
START = datetime.now()


# In[]:
files = [_ for _ in os.scandir(os.path.join(parser.parse_args().pridedb_version,'PSMs_MS1_quant')) if _.name.startswith('PSMs_MS1_quant_')]


# In[]:
combo_psms = []
for f in files:
psms = pd.read_csv(f.path, low_memory=False, usecols=['PXD','filename','scan_number','psm_id','run_id','peptidoform_id',
                                                          'precursor_mass','precursor_charge','retention_time','Peak intensity'])
    psms['file_scan'] = psms.PXD + '|' + psms.filename.apply(str) + '|' + psms.scan_number.apply(str)
    psms = psms[psms['Peak intensity']>0].copy(deep=True)
    if len(psms)!=len(set(psms.file_scan)):
        print(f.path)
        print(psms.shape)
        print(len(set(psms.file_scan)))
    combo_psms.append(psms)


# In[]:
combo_psms = pd.concat(combo_psms, ignore_index=True)


# In[]:
print('combo_psms.shape =',combo_psms.shape)
print('# unique PXD-file-scan =',len(set(combo_psms.file_scan)))
print('Any duplicates?',len(combo_psms)!=len(set(combo_psms.file_scan)))


# In[]:
combo_psms.drop(columns=['file_scan'], inplace=True)
combo_psms.to_csv(os.path.join(parser.parse_args().pridedb_version,'PSMs_MS1_quant_combined.csv.gz'), index=False, compression='gzip')


# In[]:
END = datetime.now()


# In[]:
print('Started: ', START.isoformat())
print('Finished:', END.isoformat())
