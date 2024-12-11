#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import gzip, pickle, os, re, ast, sys
from PTMmap import Fasta, PTMs_remapping
from string import ascii_uppercase
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('data_folder', metavar='data_folder', type=str, help="Path to processed folder with the data.")
parser.add_argument('date', metavar='date', type=str, help="Date in YYYY-MM-DD format.")

def ast_literal(x):
    try:
        return ast.literal_eval(x)
    except:
        return [np.nan] 

def add_seq_start(row):
    try:
        positions = np.array(row.ptm_loc)
        positions = positions + row.seq_start
        return list(positions)
    except:
        return [np.nan]


# In[ ]:
sumcounts = pd.read_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_Peptides_Abs_Counts.csv.gz")
print(sumcounts.shape)
# sumcounts.file_name = sumcounts.file_name.str.srip('.mgf.gzip')
sumcounts.reset_index(drop=True, inplace=True)
for C in ['ptm_loc','ptm_name','ptm_res','classification']:
    sumcounts[C] = sumcounts[C].apply(ast_literal)
sumcounts = sumcounts[sumcounts.ptm_loc!='[nan]'].copy(deep=True)
print(sumcounts.shape)


# In[ ]:
# sumcounts['mapped_ptms'] = sumcounts.apply(add_seq_start, axis=1)
# sumcounts.to_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_PTM_Relative_Counts_pt1_.csv.gz",
#                  index=False, compression='gzip')
xx = sumcounts.shape[0]
modcounts = []
for _,row in sumcounts.iterrows():
    print(f'{_+1:,} / {xx:,}')
    # iterator = zip(row.mapped_ptms, row.ptm_name, row.ptm_res, row.classification)
    iterator = zip(row.ptm_loc, row.ptm_name, row.ptm_res, row.classification)
    for p,m,r,c in iterator:
        if r not in list(ascii_uppercase) or c in ['ragging','semi_tryptic']:
            continue
        tmp = [
            row.file_name, row.peptidoform_id, row.psm_counts, row.is_modified, row.sequence,
            row.TOT_seq_counts, row.LeadProt, row.all_UniAcc, row.TOT_seq_counts,
            p, m, r, c
        ]
        modcounts.append(tmp)
del sumcounts

cols = [
    'file_name',
    'peptidoform_id',
    'psm_counts',
    'is_modified',
    'sequence',
    'TOT_seq_counts',
    'LeadProt',
    'Gene_name',
    'Uniprot_entry_name',
    'ptm_loc',
    'ptm_name',
    'ptm_res',
    'ptm_class',
]

modcounts = pd.DataFrame(modcounts, columns=cols)


# In[ ]:
modcounts.reset_index(drop=True, inplace=True)
modcounts.ptm_loc = modcounts.ptm_loc.apply(int)
modcounts['PTM_ID']=modcounts.apply(lambda r: f"{r.LeadProt}|{r.ptm_loc}|{r.ptm_res}|{r.ptm_name}", axis=1)
# modcounts.to_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_PTM_Relative_Counts_pt2_.csv.gz",
#                  index=False, compression='gzip')


# In[ ]:
print("...Calculating relative counts...")
relative_counts = []
for (file,ptm),df in modcounts.groupby(['file_name','PTM_ID']).__iter__():
    mod_site_peps = df.psm_counts.sum()
    tot_site_peps = df.TOT_seq_counts.sum()
    R = mod_site_peps / tot_site_peps
    relative_counts.append([file,ptm,mod_site_peps,tot_site_peps,round(R,3)])

relative_counts = pd.DataFrame(relative_counts, 
                    columns=['file_name','PTM_ID',
                             'mod_site_psms','tot_site_psms',
                             'relative_counts'])


# In[ ]:
relative_counts[['UniAcc','ptm_pos','ptm_res','ptm_name']] = relative_counts.PTM_ID.str.split('|',expand=True)
relative_counts.to_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_PTM_Relative_Counts.csv.gz",
                       index=False, compression='gzip')
