#!/usr/bin/env python3
# Test haplotypes using plink 1.07
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 10, 2023
# assoc_hap_Hdbscan v0.1
# Please check README.txt
#        check file and table names
# python3 assoc_hap_Hdbscan.py

import os

main_fp = '/mnt/wd/nsap/'

### Clusterization parameters
# eps = 0.75
# csm = 0.65
# eps_str = str(eps).replace('.', '')
# csm_str = str(csm).replace('.', '')
# min_samples = 2

# num_chr = input('Num_chr: ')

# filename = f"A_chr{num_chr}_CSM{csm_str}MS{min_samples}_EPS{eps_str}MS{min_samples}"
def assoc_hap(num_chr, filename):
  os.system(f"""{main_fp}plink-1.07-x86_64/./plink --noweb \
                        --bfile {main_fp}imp1/chr{num_chr} \
                        --hap {main_fp}Clustering/main_pipeline/hlists/{filename}_dbscan.hlist \
                        --hap-assoc \
                        --allow-no-sex \
                        --out {main_fp}Clustering/main_pipeline/hap_assoc/dbscan/{filename}_dbscan""") 
  os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/dbscan/{filename}_dbscan.log")
  os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/dbscan/{filename}_dbscan.mishap")
  try:
    os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/dbscan/{filename}_dbscan.nosex")
  except FileNotFoundError:
    print('No .nosex file found.')
  
  os.system(f"""{main_fp}plink-1.07-x86_64/./plink --noweb \
                        --bfile {main_fp}imp1/chr{num_chr} \
                        --hap {main_fp}Clustering/main_pipeline/hlists/{filename}_hdbscan.hlist \
                        --hap-assoc \
                        --allow-no-sex \
                        --out {main_fp}Clustering/main_pipeline/hap_assoc/hdbscan/{filename}_hdbscan""") 
  os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/hdbscan/{filename}_hdbscan.log")
  os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/hdbscan/{filename}_hdbscan.mishap")
  try:
    os.remove(f"{main_fp}Clustering/main_pipeline/hap_assoc/hdbscan/{filename}_hdbscan.nosex")
  except FileNotFoundError:
    print('No .nosex file found.')
    
  print(f"Haplotype restoring and testing for chromosome {num_chr} is done!")
