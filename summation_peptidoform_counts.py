#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import os, sys
from datetime import datetime
dataset, year, month, day = sys.argv[1:5]
os.mkdir(f"./{dataset}_{year}{month}{day}/tmp") 
start = datetime.now().isoformat()

# read in data
data = pd.read_csv(f"./{dataset}_{year}{month}{day}/{year}-{month}-{day}_Peptidoforms_counts_mapped.csv.gz")
data['TOT_seq_counts'] = data.psm_counts
print("data.shape =",data.shape)

# print unique peptides (those that don't need to be summed bcs they are unique)
data[~data.duplicated(['file_name','sequence'], 
                      keep=False)].to_csv(f"./{dataset}_{year}{month}{day}/tmp/{year}-{month}-{day}_Normalized_Counts_uniq.csv.gz", 
                                          index=False, compression='gzip')

data = data[data.duplicated(['file_name','sequence'], keep=False)].copy(deep=True)
data.drop(columns=['TOT_seq_counts'], inplace=True)

# lculate total psm counts per pep sequence
normalized_counts, N = [], 0
for _,df in data.groupby(['file_name','sequence']).__iter__():
    df['TOT_seq_counts'] = np.sum(df['psm_counts']) 
    normalized_counts.append(df)
    N += 1
    if N % 10000 == 0:
        print(_,f"{N:,}")
        normalized_counts = pd.concat(normalized_counts, ignore_index=True)
        normalized_counts.to_csv(f"./{dataset}_{year}{month}{day}/tmp/{year}-{month}-{day}-Normalized_Counts_{N}.csv.gz", 
                                index=False, compression='gzip')
        del normalized_counts
        normalized_counts = []

normalized_counts = pd.concat(normalized_counts, ignore_index=True)
normalized_counts.to_csv(f"./{dataset}_{year}{month}{day}/tmp/{year}-{month}-{day}-Normalized_Counts_{N}.csv.gz", 
                        index=False, compression='gzip')
del normalized_counts

print("#Peptides =",N)

final_df = [pd.read_csv(F.path) for F in os.scandir(f"./{dataset}_{year}{month}{day}/tmp/")]
# final_df.append(pd.read_csv(f"./{dataset}_{year}{month}{day}/{year}-{month}-{day}_Normalized_Counts_uniq.csv.gz"))
final_df = pd.concat(final_df, ignore_index=True)
final_df.to_csv(f"./{dataset}_{year}{month}{day}/{year}-{month}-{day}_Normalized_Counts.csv.gz", 
                index=False, compression='gzip')


print('STARTED:', start)
print('FINISHED:', datetime.now().isoformat())

