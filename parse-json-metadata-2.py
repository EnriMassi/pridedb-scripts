#!/usr/bin/env python
# coding: utf-8
import json, os, sys, re
import pandas as pd
import numpy as np

def do_command(CMD):
    print(">",CMD)
    p = os.popen(CMD).read().strip()
    print(p)
    return p


rawfiles_stats = pd.read_csv(sys.argv[1])
# output_file_name = os.path.split(sys.argv[1])[-1]
output_file_name = re.sub(r'\.csv\.gz', '_fixed.csv', sys.argv[1])

print(rawfiles_stats.tail())

def GetPSMsFirstCounts(row):
    filename = row.file_name
    pxd = row.project_id
    try:
        results_path = os.popen(f"ls -ld /public/conode5*_pride*/PRIDE_DATA/{pxd}/IONBOT_v0.11.0").read().strip().split()[-1]
        ionbot_stats = pd.read_csv(os.path.join(results_path,filename,'ionbot.stats.csv'),
                                   usecols=['#PSMs-first','#Peptides-first','#Peptidoforms-first','#PSMs-coeluting'])
        return ionbot_stats[['#PSMs-first','#Peptides-first','#Peptidoforms-first','#PSMs-coeluting']].iloc[0]
    except:
        return 'missing_results_folder!'

def GetMS2SpectraCounts(row):
    pxd = row.project_id
    try:
        spectra_path = os.popen(f"ls -ld /public/conode5*_pride*/PRIDE_DATA/{pxd}/metadata").read().strip().split()[-1]
    except:
        return 'missing_metadata_folder!'
    filename = row.file_name
    json_path = os.path.join(spectra_path, re.sub(r'\.mgf\.gzip','-metadata.json',filename))

    try:
        with open(json_path) as user_file:
            parsed_json = json.load(user_file)
    except:
        return 'missing_metadata_file!'
    
    for entry in parsed_json['MsData']:
        if entry["name"] == "Number of MS2 spectra":
            return int(entry['value'])

def isHCD(row):
    pxd = row.project_id
    try:
        spectra_path = os.popen(f"ls -ld /public/conode5*_pride*/PRIDE_DATA/{pxd}/metadata").read().strip().split()[-1]
    except:
        return 'missing_metadata_folder!'
    filename = row.file_name
    json_path = os.path.join(spectra_path, re.sub(r'\.mgf\.gzip','-metadata.json',filename))

    try:
        with open(json_path) as user_file:
            parsed_json = json.load(user_file)
    except:
        return 'missing_metadata_file!'
    
    for entry in parsed_json['ScanSettings']:
        if entry["value"] == "HCD":
            return True
    return False

def isCID(row):
    pxd = row.project_id
    try:
        spectra_path = os.popen(f"ls -ld /public/conode5*_pride*/PRIDE_DATA/{pxd}/metadata").read().strip().split()[-1]
    except:
        return 'missing_metadata_folder!'
    filename = row.file_name
    json_path = os.path.join(spectra_path, re.sub(r'\.mgf\.gzip','-metadata.json',filename))

    try:
        with open(json_path) as user_file:
            parsed_json = json.load(user_file)
    except:
        return 'missing_metadata_file!'
    
    for entry in parsed_json['ScanSettings']:
        if entry["value"] == "CID":
            return True
    return False
    
rawfiles_stats[['#PSMs-first','#Peptides-first','#Peptidoforms-first','#PSMs-coeluting']] = rawfiles_stats.apply(GetPSMsFirstCounts, axis=1, result_type='expand')
rawfiles_stats['Total_msms_in_file'] = rawfiles_stats.apply(GetMS2SpectraCounts, axis=1)
# rawfiles_stats['isHCD'] = rawfiles_stats.apply(isHCD, axis=1)
# rawfiles_stats['isCID'] = rawfiles_stats.apply(isCID, axis=1)
rawfiles_stats.to_csv(output_file_name, index=False)
print(output_file_name)