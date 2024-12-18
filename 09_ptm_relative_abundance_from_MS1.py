#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import ast, os
from datetime import datetime
from string import ascii_uppercase

def ast_literal(x):
    try:
        return ast.literal_eval(x)
    except:
        return [np.nan] 

def zip_modifications(row):
    iterator = zip(row.ptm_loc, row.ptm_name, row.ptm_res, row.classification)
    return [f'{row.LeadProt}|{p}|{r}|{m}' for p,m,r,c in iterator if r in list(ascii_uppercase) and c not in ['ragging','semi_tryptic']]

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('data_folder', metavar='data_folder', type=str, help="Path to processed folder with the data.")
parser.add_argument('date', metavar='date', type=str, help="Date in YYYY-MM-DD format.")

START = datetime.now()


# In[4]:
print('> Reading peptidoforms...')
peptidoforms = pd.read_csv(os.path.join(parser.parse_args().data_folder,
                                        f'{parser.parse_args().date}_Unique_peptidoforms.csv.gz'),
                           low_memory=False)
for C in ['ptm_loc','ptm_name','ptm_res','classification']: #
    peptidoforms[C] = peptidoforms[C].apply(ast_literal)
peptidoforms.ptm_loc = peptidoforms.ptm_loc.apply(lambda x: [np.nan] if x==[np.nan] else [str(int(_)) for _ in x])

peptidoforms['modifications'] = peptidoforms.apply(zip_modifications, axis=1)
peptidoforms.drop(columns=['ptm_name','ptm_loc','ptm_res','classification'], inplace=True)

peptidoforms = peptidoforms.explode('modifications')


# In[6]:
print('> Reading peak intensities...')
data = pd.read_csv(os.path.join(parser.parse_args().data_folder,
                                f'{parser.parse_args().date}_PSMs_MS1_median_peak_intensities.csv.gz'), 
                   low_memory=False)


# In[7]:
print('> merging data frames...')
merged_data = data.merge(peptidoforms, on='peptidoform_id')
print(merged_data.shape)
del data, peptidoforms


# In[8]:
# sanity check
print(merged_data.drop_duplicates(['PXD','filename','peptidoform_id','modifications']).shape)


# In[9]:
print('> Normalizing...')
normalized = []
for _,df in merged_data.groupby(['PXD','filename']):
    df.log_peak_intensity -= np.min(df.log_peak_intensity)
    df.log_peak_intensity /= np.max(df.log_peak_intensity)
    normalized.append(df)
del df
normalized = pd.concat(normalized, ignore_index=True) 
normalized.shape


# In[11]:
normalized[['LeadProt2','ptm_loc','ptm_res','ptm_name']] = normalized.modifications.str.split('|', expand=True)
normalized


# In[12]:
normalized['isContaminant'] = normalized.LeadProt.apply(lambda x: 'CONTAMINANT' in x.upper())
normalized = normalized[~normalized.isContaminant].copy(deep=True)
normalized.shape


# In[15]:
print('> Calculating medians...')
medians_per_ptm = normalized.groupby(['PXD','filename','LeadProt','ptm_name','ptm_loc','ptm_res'])[['log_peak_intensity']].median()
medians_per_ptm.columns = ['median_normal_intensity']
medians_per_ptm.reset_index(inplace=True)
print('medians_per_ptms.shape =',medians_per_ptm.shape)


# In[16]:
print('> printing to', os.path.join(parser.parse_args().data_folder, f'{parser.parse_args().date}_PTMs_intensities.csv.gz'))
medians_per_ptm.to_csv(os.path.join(parser.parse_args().data_folder, f'{parser.parse_args().date}_PTMs_intensities.csv.gz'), 
                       index=False, compression='gzip')

# In[]:
END = datetime.now()
print('Started: ', START.isoformat())
print('Finished:', END.isoformat())