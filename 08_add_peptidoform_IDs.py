#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import ast, os
from datetime import datetime
START = datetime.now()

def ast_literal(x):
    try:
        return ast.literal_eval(x)
    except:
        return [np.nan] 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('peptidoforms_file', metavar='peptidoforms_file', type=str, help="Path to file with peptidoforms data.")
parser.add_argument('other_file', metavar='other_file', type=str, help="Path to file with intensitites or psms data (which peptidoform data will be added to).")
parser.add_argument('output', metavar='output', type=str, help="name of the output file (without extension!). This script can do 2 types of merge.")


# In[5]:
peptidoforms = pd.read_csv(parser.parse_args().peptidoforms_file, low_memory=False, 
                            # usecols=[1,4,5,6,7,8,9,10,11,12]
                        )
for C in ['ptm_loc','ptm_name','ptm_res','classification']: #
    peptidoforms[C] = peptidoforms[C].apply(ast_literal)
peptidoforms.ptm_loc = peptidoforms.ptm_loc.apply(lambda x: [np.nan] if x==[np.nan] else [str(int(_)) for _ in x])
for C in ['ptm_loc','ptm_name','ptm_res','classification']: #
    peptidoforms[C] = peptidoforms[C].apply(lambda x: np.nan if x==[np.nan] else '|'.join(x))
# peptidoforms.drop_duplicates(inplace=True)
print(f'Length of peptidoforms df = {peptidoforms.shape[0]:,}')
print(f'# unique peptidoform IDs = {len(set(peptidoforms.peptidoform_id)):,}')
print('-------')


# In[5]:
data = pd.read_csv(parser.parse_args().other_file, low_memory=False)
merged_data = peptidoforms.merge(data, on='peptidoform_id')

# sanity check
try:
    print(merged_data.drop_duplicates(['PXD','filename','peptidoform_id']).shape, merged_data.shape)
except:
    print(f'Length of peptidoforms df = {merged_data.shape[0]:,}')
    print(f'# unique peptidoform IDs = {len(set(merged_data.psm_id)):,}')

# merged_data.to_csv(parser.parse_args().output+'.csv.gz', index=False, compression='gzip')


# In[]:
END = datetime.now()
print('Started: ', START.isoformat())
print('Finished:', END.isoformat())