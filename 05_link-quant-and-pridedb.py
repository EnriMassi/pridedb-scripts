#!/usr/bin/env python
# coding: utf-8
import os, sys
import pandas as pd
import numpy as np
from datetime import datetime
START = datetime.now()

def do_command(CMD):
    print(">",CMD)
    p = os.popen(CMD).read().strip()
    print(p)

def read_QuantifiedPeaks_file(PXD_):
    CMD = f'ls -ld /public/conode5*_pride*/PRIDE_DATA/{PXD_}/QuantifiedPeaks_v2.tsv'
    p = os.popen(CMD).read().strip()
    quant_path = p.split()[-1]
    x = pd.read_csv(quant_path, sep='\t', low_memory=False)
    x.dropna(axis=1, inplace=True, how='all')
    x.rename(columns={'Protein Group':'spectrum_title', 'File Name':'filename'}, inplace=True)
    x.spectrum_title = x.spectrum_title.str.split(';')
    x = x.explode('spectrum_title')
    x['scan_number'] = x.spectrum_title.apply(lambda x: x.split('scan=')[-1])
    x.drop(columns=['Base Sequence','Full Sequence'], inplace=True)
    x['PXD'] = PXD_
    x = x.set_index(['PXD','filename','scan_number']).reset_index()
    x = x[x['Peak RT Apex']!='-'].copy(deep=True)

    # set columns types and round floats
    for C in ['Peptide Monoisotopic Mass','MS2 Retention Time','Theoretical MZ', 
              'Peak intensity','Peak RT Start','Peak RT Apex','Peak RT End','Peak MZ','Peak Charge',
              'Peak Split Valley RT','Peak Apex Mass Error (ppm)']:
        x[C] = x[C].apply(float).round(4)
    for C in ['Peak Charge']:
        x[C] = x[C].apply(int)

    return x

# In[]
try:
    psms_file = sys.argv[1]
except:
    psms_file = '/home/enrico/FWO_project/MoDPA/v0113/NonDuplicated_PSMs_v0113.csv.gz'

file_fld = os.path.split(psms_file)[0]
os.makedirs(os.path.join(file_fld,'PSMs_MS1_quant'))

psms = pd.read_csv(psms_file, low_memory=False)
psms.scan_number = psms.scan_number.apply(int)
psms.PXD         = psms.PXD.apply(str)
psms.filename    = psms.filename.apply(str)
print(f'#PSMs = {len(psms):,}','\n')

# # In[]
pxds = list(set(psms.PXD))
pxds.sort()
skipped_pxds = []
for i,PXD in enumerate(pxds):
    filepath = os.path.join(file_fld,'PSMs_MS1_quant',f'PSMs_MS1_quant_{PXD}.csv.gz')
    if os.path.exists(filepath):
        continue
    print(PXD, '--', i+1, '/', len(pxds))
    try:
        quant = read_QuantifiedPeaks_file(PXD)   
    except:
        skipped_pxds.append(PXD)
        continue
    quant.scan_number   = quant.scan_number.apply(int)
    quant.PXD           = quant.PXD.apply(str)
    quant.filename      = quant.filename.apply(str)
    psms.precursor_mass = psms.precursor_mass.apply(float).round(4)
    psms.retention_time = psms.retention_time.apply(float).round(4)
    combo = psms.merge(quant, #[['PXD','filename','scan_number','Peak intensity']],
                       on=['PXD','filename','scan_number'])
    print(len(combo))
    print(len(set(combo['psm_id'])))
    print(len(combo.drop_duplicates()))
    if len(combo) > 0:
        combo.to_csv(filepath, index=False, compression='gzip')
    del quant, combo
print('-----------')

with open(os.path.join(file_fld,'PSMs_MS1_quant','skipped-projects.txt'), 'w') as OUT:
    for j in skipped_pxds:
        OUT.write(j)
        OUT.write('\n')

END = datetime.now()
print('\n')
print(START)
print(END)