#!/usr/bin/env python3
# Consolidate assoc.hap files, sort entries by p-value, wirte csv outfile
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 19, 2023
# CHR_START_check.py
# Please check README.txt
#        check prefixes
# python3 CHR_START_check.py

import pandas as pd
import os
import re

# folder path
main_fp = '/mnt/wd/nsap/Clustering/main_pipeline/hap_assoc/csvs/'
# DBSCAN
# Specify the prefix for omnibus_not_filtered_
file_prefix = 'omnibus_not_filtered_A_chr'
matching_dbscan = [
    file_name
    for file_name in os.listdir(main_fp + 'dbscan/')
    if file_name.startswith(file_prefix)
]
# for not filtered omnibus files
counter = 0
for dbscan_name in matching_dbscan:
  # reading file
  nOmni_pd = pd.read_csv(main_fp + 'dbscan/' + dbscan_name)
  # forming a dictionary with unique CHR:START keys and appending END to a list under that key
  nOmni_counter = {}
  # 'LOCUS' contains CHR:START-END
  for row in nOmni_pd['LOCUS']:
    # key will be CHR:START
    key = str(row).split('-')[0]
    if not key in nOmni_counter:
            nOmni_counter[key] = []
    nOmni_counter[key].append(str(row).split('-')[1])
  # print(len(nOmni_counter))
  counter += len(nOmni_counter)
print('Totall number of cluster in \'omnibus_not_filtered_\' DBSCAN files is: ', counter)

# Specify the prefix and postfix for omnibus
file_prefix = 'A_chr'
file_postfix = 'omnibus.csv'
matching_dbscan = [
    file_name
    for file_name in os.listdir(main_fp + 'dbscan/')
    if file_name.startswith(file_prefix) and file_name.endswith(file_postfix)
]
# for omnibus files
counter = 0
for dbscan_name in matching_dbscan:
  # reading file
  omni_pd = pd.read_csv(main_fp + 'dbscan/' + dbscan_name)
  # forming a dictionary with unique CHR:START keys and appending END to a list under that key
  omni_counter = {}
  # 'LOCUS' contains CHR:START-END
  for row in omni_pd['LOCUS']:
    # key will be CHR:START
    key = str(row).split('-')[0]
    if not key in omni_counter:
            omni_counter[key] = []
    omni_counter[key].append(str(row).split('-')[1])
  # print(len(omni_counter))
  counter += len(omni_counter)
print('Totall number of cluster in \'_omnibus_\' DBSCAN files is: ', counter)

# HDBSCAN
# Specify the prefix for omnibus_not_filtered_
file_prefix = 'omnibus_not_filtered_A_chr'
matching_hdbscan = [
    file_name
    for file_name in os.listdir(main_fp + 'hdbscan/')
    if file_name.startswith(file_prefix)
]
# for not filtered omnibus files
counter = 0
for hdbscan_name in matching_hdbscan:
  # reading file
  nOmni_pd = pd.read_csv(main_fp + 'hdbscan/' + hdbscan_name)
  # forming a dictionary with unique CHR:START keys and appending END to a list under that key
  nOmni_counter = {}
  # 'LOCUS' contains CHR:START-END
  for row in nOmni_pd['LOCUS']:
    # key will be CHR:START
    key = str(row).split('-')[0]
    if not key in nOmni_counter:
            nOmni_counter[key] = []
    nOmni_counter[key].append(str(row).split('-')[1])
  # print(len(nOmni_counter))
  counter += len(nOmni_counter)
print('Totall number of cluster in \'omnibus_not_filtered_\' HDBSCAN files is: ', counter)

# Specify the prefix and postfix for omnibus
file_prefix = 'A_chr'
file_postfix = 'omnibus.csv'
matching_hdbscan = [
    file_name
    for file_name in os.listdir(main_fp + 'hdbscan/')
    if file_name.startswith(file_prefix) and file_name.endswith(file_postfix)
]
# for omnibus files
counter = 0
for hdbscan_name in matching_hdbscan:
  # reading file
  omni_pd = pd.read_csv(main_fp + 'hdbscan/' + hdbscan_name)
  # forming a dictionary with unique CHR:START keys and appending END to a list under that key
  omni_counter = {}
  # 'LOCUS' contains CHR:START-END
  for row in omni_pd['LOCUS']:
    # key will be CHR:START
    key = str(row).split('-')[0]
    if not key in omni_counter:
            omni_counter[key] = []
    omni_counter[key].append(str(row).split('-')[1])
  # print(len(omni_counter))
  counter += len(omni_counter)
print('Totall number of cluster in \'_omnibus_\' HDBSCAN files is: ', counter)
