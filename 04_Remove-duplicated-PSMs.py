#!/usr/bin/env python
# coding: utf-8
import os, sys
import pandas as pd
import numpy as np
from datetime import datetime
START = datetime.now()
filepath = sys.argv[1]


# In[3]:
def read_PSMs_file(path_to_file):
    df = pd.read_csv(path_to_file, low_memory=False)
    df.usi = df.usi.str.split(':')
    df['scan_number'] = df.usi.apply(lambda x: x[4])
    df['PXD']         = df.usi.apply(lambda x: x[1])
    df['filename']    = df.usi.apply(lambda x: x[2].strip('.mgf'))
    df.drop(columns=['usi'], inplace=True)
    df.set_index(['PXD','filename','scan_number'], inplace=True)
    df.reset_index(inplace=True)
    df = df.reset_index(drop=True)
    return df


# In[4]:
psms = read_PSMs_file(filepath)
print(f'{len(psms):,}')
psms2 = psms[~psms.duplicated(['PXD','filename','scan_number'],keep=False)].copy(deep=True)
print(f'{len(psms2):,}')

# In[5]:
folder,file = os.path.split(filepath)
psms2.to_csv(os.path.join(folder,file.replace('_PSMs_','_PSMs_unique_')), index=False)

# In[]:
psms2 = psms[psms.duplicated(['PXD','filename','scan_number'],keep=False)].copy(deep=True)
psms2.to_csv(os.path.join(folder,file.replace('_PSMs_','_PSMs_duplicated_')))
print(f'{len(psms2):,}')

# In[]:
duplications_counts = psms[['PXD','filename','scan_number']].value_counts().reset_index()
duplications_counts.columns = ['PXD','filename','scan_number','dup_counts']
X = duplications_counts.dup_counts.value_counts(normalize=True).round(2)
print(X)