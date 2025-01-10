#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import os, sys
from datetime import datetime
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('data_folder', metavar='data_folder', type=str, help="Path to processed folder with the data.")
parser.add_argument('date', metavar='date', type=str, help="Date in YYYY-MM-DD format.")

START = datetime.now()
# os.mkdir(os.path.join(parser.parse_args().data_folder,'tmp')) 

# read in data
data = pd.read_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_Peptidoforms_counts_mapped.csv.gz")

# In[new version -- test]
print("data.shape =",data.shape)
tmp = data.groupby(['file_name','sequence'])[['psm_counts']].sum()
tmp.columns = ['TOT_seq_counts']
tmp.reset_index(inplace=True)
tmp = data.merge(tmp, on=['file_name','sequence'])
print("#Peptides =", len(set(tmp.sequence)))
tmp.to_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_Peptides_Abs_Counts.csv.gz", 
           index=False, compression='gzip')

# In[old version]
# data['TOT_seq_counts'] = data.psm_counts
# print("data.shape =",data.shape)

# # print unique peptides (those that don't need to be summed bcs they are unique)
# data[~data.duplicated(['file_name','sequence'], 
#                       keep=False)].to_csv(f"./{parser.parse_args().data_folder}/tmp/{parser.parse_args().date}_Peptides_Abs_Counts_uniq.csv.gz", 
#                                           index=False, compression='gzip')

# data = data[data.duplicated(['file_name','sequence'], keep=False)].copy(deep=True)
# data.drop(columns=['TOT_seq_counts'], inplace=True)

# # lculate total psm counts per pep sequence
# normalized_counts, N = [], 0
# for _,df in data.groupby(['file_name','sequence']).__iter__():
#     df['TOT_seq_counts'] = np.sum(df['psm_counts']) 
#     normalized_counts.append(df)
#     N += 1
#     if N % 10000 == 0:
#         print(_,f"{N:,}")
#         normalized_counts = pd.concat(normalized_counts, ignore_index=True)
#         normalized_counts.to_csv(f"./{parser.parse_args().data_folder}/tmp/{parser.parse_args().date}_Peptides_Abs_Counts_{N}.csv.gz", 
#                                 index=False, compression='gzip')
#         del normalized_counts
#         normalized_counts = []

# normalized_counts = pd.concat(normalized_counts, ignore_index=True)
# normalized_counts.to_csv(f"./{parser.parse_args().data_folder}/tmp/{parser.parse_args().date}_Peptides_Abs_Counts_{N}.csv.gz", 
#                         index=False, compression='gzip')
# del normalized_counts

# print("#Peptides =",N)

# final_df = [pd.read_csv(F.path) for F in os.scandir(f"./{parser.parse_args().data_folder}/tmp/")]
# final_df = pd.concat(final_df, ignore_index=True)
# final_df.to_csv(f"./{parser.parse_args().data_folder}/{parser.parse_args().date}_Peptides_Abs_Counts.csv.gz", 
#                 index=False, compression='gzip')



# In[]:
END = datetime.now()
print('Started: ', START.isoformat())
print('Finished:', END.isoformat())