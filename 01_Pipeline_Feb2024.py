#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
from PTMmap import Fasta
from datetime import date
TODAY = date.today().isoformat() 
print('The date is:',TODAY)


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('peptidoform_ids', metavar='peptidoform_ids', type=str, 
                    help="Path to 'Peptidoforms' file.")
parser.add_argument('peptidoform_counts', metavar='peptidoform_counts', type=str, 
                    help="Path to 'peptidoform_counts' file.")
parser.add_argument('peptides_mappings', metavar='peptides_mappings', type=str, 
                    help="Path to 'peptides_mappings' file.")
parser.add_argument('fasta', metavar='fasta', type=str, 
                    help='Path to FASTA file.')

fasta = {}
Fasta.getFasta(parser.parse_args().fasta, fasta)
Fasta.addClassification(fasta)

# In[]: to sort canonical prots, isoforms, ORFs, etc.
def getClass(_):
    try:
        return fasta[_]['Class']
    except:
        return 'zzz'
getClassVect = np.vectorize(getClass)


# In[2]:
def read_peptide_to_protein_mappings(pep_dict_path, pep_set):
    # read ionbot pep dict file
    pepdict = pd.read_csv(pep_dict_path)
    print("pepdict size:",pepdict.shape)
    print(f"#Unique pep sequences in searchDB = {len(set(pepdict.peptide)):,}")
    pepdict = pepdict[pepdict.peptide.isin(pep_set)].copy(deep=True)
    print(f"#Unique pep sequences found by Ionbot = {len(set(pepdict.peptide)):,}")

    # rearrange part 1
    pepdict.proteins = pepdict.proteins.str.split('\|\|')
    pepdict = pepdict.explode('proteins')
    pepdict.proteins = pepdict.proteins.str.split('\(\(|\)\)')
    pepdict['UniAcc'] = pepdict.proteins.apply(lambda x: x[3].replace('|','.'))
    pepdict['entry']  = pepdict.proteins.apply(lambda x: x[0])
    pepdict[['start','end']]=pepdict.proteins.apply(lambda x: x[1]).str.split('-', expand=True).astype(int)

    # rearrange part 2
    maps_ambig = []
    for pep_seq,df in pepdict.groupby('peptide'):
        df.sort_values('UniAcc', inplace=True, key=getClassVect)
        maps_ambig.append([pep_seq, list(df.start), list(df.UniAcc), list(df.entry)])
    maps_ambig = pd.DataFrame(maps_ambig, columns=['sequence','all_pep_start','all_UniAcc','all_Entries'])
    maps_ambig['ambiguous_map'] = maps_ambig.all_UniAcc.apply(lambda x: len(x)>1)
    return maps_ambig


# In[3]:
def unique_counts_per_UniAcc(prot_list, counts):
    return [counts[_] if _ in list(counts.keys()) else 0 for _ in prot_list]


# In[4]:
def get_leading(row):
    counts = row['unique_counts']
    accs = row['all_UniAcc']
    entries = row['all_Entries']
    starts = row['all_pep_start']
    leading_prot = accs[counts.index(max(counts))]
    leading_entry = entries[counts.index(max(counts))]
    leading_start = starts[counts.index(max(counts))]
    return leading_prot, leading_entry, leading_start


# In[5]:
def get_ambiguous_peptides(pep_dict_path, pep_set):
    maps_ambig_partial = read_peptide_to_protein_mappings(pep_dict_path, pep_set)
    # maps_ambig_partial = maps_ambig_partial[maps_ambig_partial.sequence.isin(pep_set)]
    # print(f"#Unique pep sequences found by Ionbot = {len(set(maps_ambig_partial.sequence)):,}")
    
    unique_mappings = maps_ambig_partial[~maps_ambig_partial.ambiguous_map].copy(deep=True)
    unique_mappings['unique_protein'] = unique_mappings.all_UniAcc.apply(lambda x: x[0])
    # unique_mappings.sort_values('unique_protein')
    
    unique_counts = unique_mappings.groupby('unique_protein').count()[['sequence']]
    unique_counts = unique_counts.to_dict()['sequence']
    # unique_counts
    
    maps_ambig_partial['unique_counts'] = maps_ambig_partial.all_UniAcc.apply(lambda x: unique_counts_per_UniAcc(x,unique_counts))
    
    maps_ambig_partial[['LeadProt','LeadEntry','pep_start']] = maps_ambig_partial.apply(get_leading, axis=1, 
                                                                                        result_type='expand')
    maps_ambig_partial.all_UniAcc = maps_ambig_partial.all_UniAcc.apply(lambda x: '||'.join(x))
    maps_ambig_partial.sort_values('ambiguous_map', inplace=True)
    return maps_ambig_partial[['sequence','pep_start','LeadProt','LeadEntry','all_UniAcc']]


# In[6]:
def map_ionbot_IDs(IDs_path, pep_dict_path):
    ids = pd.read_csv(IDs_path)
    
    maps_ambig = get_ambiguous_peptides(pep_dict_path,set(ids.sequence))
    
    print("IDs table size:", ids.shape)
    mapped_ids = ids.merge(maps_ambig, on='sequence')
    print("Mapped_IDs table size:",mapped_ids.shape)
    mapped_ids.ptm_loc = mapped_ids.ptm_loc + mapped_ids.pep_start - 1
    mapped_ids.drop(columns=['pep_start'], inplace=True)
    
    return mapped_ids


# In[7]:
def group_IDs_into_peptidoforms(mapped_ids):
    mapped_ids_grouped = []
    iterator = mapped_ids.groupby(['peptidoform_id','peptide_id','is_modified',
                                   'sequence','LeadProt','LeadEntry','all_UniAcc'])
    for (peptidoform_id,peptide_id,is_modified,seq,leadprot,leadentry,allprots),df in iterator.__iter__():
        mapped_ids_grouped.append([
            peptidoform_id,
            peptide_id,
            is_modified=='t',
            seq,
            list(df.ptm_name),
            list(df.ptm_loc),
            list(df.ptm_res),
            list(df.classification),
            leadprot,
            leadentry,
            allprots
        ])
    return pd.DataFrame(mapped_ids_grouped, columns=mapped_ids.columns)


# In[8]:
def add_psm_counts(mapped_peptidoforms, psm_counts_path):
    counts  = pd.read_csv(psm_counts_path, usecols=['file_name','peptidoform_id','psm_counts'])
    print(counts.shape)
    print(mapped_peptidoforms.shape)
    counts = counts.merge(mapped_peptidoforms, on='peptidoform_id')
    print(counts.shape)
    return counts


# In[9]:
def psm_counts_per_PTM(mapped_peptidoforms_counts):
    modcounts = []
    for _,row in mapped_peptidoforms_counts.iterrows():
        iterator = zip(row.ptm_loc, row.ptm_name, row.ptm_res, row.classification)
        for p,m,r,c in iterator:
            if c=='ragging':
                continue
            tmp = [
                row.file_name, row.peptidoform_id, row.psm_counts, row.sequence, row.is_modified,  
                row.LeadProt, row.LeadEntry, row.all_UniAcc,
                p, m, r, c
            ]
            modcounts.append(tmp)
    
    cols = [
        'file_name',
        'peptidoform_id',
        'psm_counts',
        'sequence',
        'is_modified',
        'LeadProt',
        'LeadEntry',
        'all_UniAcc',
        'ptm_loc',
        'ptm_name',
        'ptm_res',
        'ptm_class'
    ]
    modcounts = pd.DataFrame(modcounts, columns=cols)
    
    modcounts = modcounts.groupby(['file_name','LeadProt','LeadEntry',
                                   'ptm_loc','ptm_name','ptm_res','ptm_class']).sum()[['psm_counts']]
    modcounts.columns = ['total_counts']
    print('\n',modcounts.describe(),'\n')
    return modcounts.reset_index()



# In[10]:
species = 'human'
mapped_ids = map_ionbot_IDs(parser.parse_args().peptidoform_ids, parser.parse_args().peptides_mappings)
peptidoforms = group_IDs_into_peptidoforms(mapped_ids)
peptidoforms.to_csv(f'{TODAY}_Unique_peptidoforms.csv.gz', compression='gzip',index=False)
mapped_peptidoforms_counts = add_psm_counts(peptidoforms, parser.parse_args().peptidoform_counts)


print(f"#peptidoforms = {len(set(peptidoforms.peptidoform_id)):,}")
print("Unique peptidoforms:", peptidoforms.shape)

mod, unmod = peptidoforms.is_modified.value_counts()
print('_')
print(peptidoforms.is_modified.value_counts())
print(f"% Unmodified peptides = {unmod/(unmod+mod):.1%}",'\n')

print("Mapped peptidoforms with counts:", mapped_peptidoforms_counts.shape)
mapped_peptidoforms_counts.to_csv(f'{TODAY}_Peptidoforms_counts_mapped.csv.gz', compression='gzip',index=False)


PTM_counts = psm_counts_per_PTM(mapped_peptidoforms_counts)
PTM_counts.to_csv(f'{TODAY}_PTM_sites_counts.csv.gz', compression='gzip',index=False)


PTM_counts['PTM_ID'] = PTM_counts.LeadProt.apply(str)+'|'+PTM_counts.ptm_res.apply(str)+'|'+PTM_counts.ptm_loc.apply(int).apply(str)+'|'+PTM_counts.ptm_name.apply(str)
PTM_counts.drop(columns=['file_name','LeadProt','ptm_loc','ptm_name','ptm_res','total_counts'], inplace=True)
# PTM_counts.drop_duplicates('PTM_ID', inplace=True)
PTM_counts.drop_duplicates(inplace=True)
print('#Unique PTMs =',PTM_counts.shape)

PTM_counts.to_csv(f'{TODAY}_ALL_PTMs.csv.gz', compression='gzip', index=False)