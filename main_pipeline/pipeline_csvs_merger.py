#!/usr/bin/env python3
# Merge csvs for all chromosomes. For files with a prefix
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 18, 2023
# Please check README.txt
# python3 pipeline_csvs_merger.py

import pandas as pd
import os
import re

# Specify the directory where the files are located
main_fp = '/mnt/wd/nsap/Clustering/main_pipeline/hap_assoc/csvs/'  
# Specify the prefix
file_prefix = 'omnibus_not_filtered_A_chr'

# DBSCAN
matching_dbscan = [
    file_name
    for file_name in os.listdir(main_fp + 'dbscan/')
    if file_name.startswith(file_prefix)
]
pd_df = pd.DataFrame({'LOCUS':[]
                     ,'HAPLOTYPE':[]
                     ,'F_A':[]
                     ,'F_U':[]
                     ,'CHISQ':[]
                     ,'DF':[]
                     ,'P':[]
                     ,'SNPS':[]
                     ,'CHR':[]
                     ,'START':[]
                     ,'END':[]
                     ,'LENGTH':[]
                     ,'SIZE':[]})
for dbscan_name in matching_dbscan:
  num_chr = re.findall("\\D*(\\d*)", dbscan_name.split('_')[4])[0]
  print(f"Working with file for {num_chr} chromosome...")
  pd_df1 = pd.read_csv(dbscan_name, lineterminator='\n', delim_whitespace=True)


# # HDBSCAN  
# matching_hdbscan = [
#     file_name
#     for file_name in os.listdir(main_fp + 'hdbscan/')
#     if file_name.startswith(file_prefix)
# ]
# for hdbscan_name in matching_hdbscan:
#   print(hdbscan_name)
