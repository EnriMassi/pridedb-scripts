#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import numpy as np
import os
from datetime import datetime
START = datetime.now()

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('pridedb_version', metavar='pridedb_version', type=str, help="Path to raw folder with the data.")
parser.add_argument('data_folder', metavar='data_folder', type=str, help="Path to processed folder with the data.")
parser.add_argument('date', metavar='date', type=str, help="Date in YYYY-MM-DD format.")


# In[3]: calculate median intensities
psms = pd.read_csv(os.path.join(parser.parse_args().pridedb_version,'PSMs_MS1_quant_combined.csv.gz'), 
                   low_memory=False, 
                   usecols=['PXD','filename','scan_number','psm_id','run_id','peptidoform_id',
                            'precursor_mass','precursor_charge','retention_time','Peak intensity'])
psms.rename(columns={'Peak intensity':'peak_intensity'}, inplace=True)
print(psms.shape)
psms = psms[psms['peak_intensity']>0].copy(deep=True)
print(psms.shape)

median_intensities = psms.groupby(['PXD','filename','peptidoform_id'])[['peak_intensity']].median()
median_intensities['log_peak_intensity'] = median_intensities.peak_intensity.apply(np.log10)
median_intensities.reset_index(inplace=True)

median_intensities.to_csv(os.path.join(parser.parse_args().data_folder,
                                       f'{parser.parse_args().date}_PSMs_MS1_median_peak_intensities.csv.gz'), 
                          index=False, compression='gzip')


# In[]:
END = datetime.now()
print('Started: ', START.isoformat())
print('Finished:', END.isoformat())
