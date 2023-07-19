#!/usr/bin/env python3
# Consolidate assoc.hap files, sort entries by p-value, wirte csv outfile
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 11, 2023
# merger_assoc_hap v0.1
# Please check README.txt
#        check file and table names
# python3 merger_assoc_hap.py -n 1 -m all

import pandas as pd
import argparse

main_fp = '/mnt/wd/nsap/Clustering/main_pipeline/'

### Clusterization parameters
eps = 0.75
csm = 0.65
eps_str = str(eps).replace('.', '')
csm_str = str(csm).replace('.', '')
min_samples = 2

def dataframe_merger(num_chr_s, num_chr_e, meth):
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
    for num_chr in range(num_chr_s, num_chr_e + 1):
      print(f"Working with file for {num_chr} chromosome...")
      ### filename
      filename = f"A_chr{num_chr}_CSM{csm_str}MS{min_samples}_EPS{eps_str}MS{min_samples}"
      
      with open(f"{main_fp}hap_assoc/{meth}/{filename}_{meth}.assoc.hap", 'r') as file:
          print(f"\rReading {num_chr} file...")
          pd_df1 = pd.read_csv(file, lineterminator='\n', delim_whitespace=True)
          ## Split LOCUS and SNPS columns into CHR, START, END, LENGTH, SIZE columns
          CHR_list, START_list, END_list, LENGTH_list, SIZE_list = [],[],[],[],[]
          for index, row in pd_df1.iterrows():
            locus_str = str(row['LOCUS']).split(':')
            CHR_list.append(locus_str[0])
            start_end = str(locus_str[1]).split('-')
            START_list.append(start_end[0])
            END_list.append(start_end[1])
            LENGTH_list.append(str(int(start_end[1]) - int(start_end[0])))
            snp_str = str(row['SNPS']).split('|')
            SIZE_list.append(len(snp_str))
          pd_df1['CHR'] = CHR_list
          pd_df1['START'] = START_list
          pd_df1['END'] = END_list
          pd_df1['LENGTH'] = LENGTH_list
          pd_df1['SIZE'] = SIZE_list
          pd_df = pd.concat([pd_df, pd_df1], sort=True)
      ### sort dataframe by P-value
      pd_df.sort_values(by=['P'], inplace=True)
      print(pd_df)
      ### creating and writing to the outfile
      ### OMNIBUS filtered
      omnibus_df = pd_df.where(pd_df['HAPLOTYPE'] == 'OMNIBUS').dropna(subset=['LOCUS'])
      with open(f"{main_fp}hap_assoc/csvs/{meth}/{filename}_{meth}_omnibus.csv", 'w') as omnibus_excel:
          omnibus_df.to_csv(omnibus_excel, lineterminator='\n', index=False)
          print('Created omnibus outfile...')
      non_omnibus_df = pd_df[pd_df.HAPLOTYPE != 'OMNIBUS']
      with open(f"{main_fp}hap_assoc/csvs/{meth}/{filename}_{meth}_nonOmnibus.csv", 'w') as non_omnibus_excel:
          non_omnibus_df.to_csv(non_omnibus_excel, lineterminator='\n', index=False)
          print('Created non omnibus outfile...')
      ### OMNIBUS not filtered
      with open(f"{main_fp}hap_assoc/csvs/{meth}/omnibus_not_filtered_{filename}_{meth}.csv", 'w') as csv_file:
          pd_df.to_csv(csv_file, lineterminator='\n', index=False)
          print('Created not filtered outfile...')
      print('Done! DataFrame is merged.')

# Parse command line arguments
parser = argparse.ArgumentParser(description='Consolidate assoc.hap files, sort entries by p-value, wirte csv outfile')
parser.add_argument('-n', '--num_chr', 
                    help = 'number of a chromosome either sole or separated with a \'-\'', 
                    required = True)
parser.add_argument('-m', '--method', 
                    help = 'possible entries are: hdbscan, dbscan, all', 
                    required = True)
# TODO search in the folder by keywords: e.g. param_*
# parser.add_argument('-i', '--input', 
#                     help = 'possible entries are: hdbscan, dbscan, all', 
#                     required = True)
args = parser.parse_args()

# Check if -n argument is sole or range
try:
  num_chr_s = int(args.num_chr.split('-')[0])
  num_chr_e = int(args.num_chr.split('-')[1])
except IndexError:
  num_chr_s = int(args.num_chr)
  num_chr_e = int(num_chr_s)

methods = args.method
if methods == 'all':
  methods = ['dbscan', 'hdbscan']
for meth in methods:
  dataframe_merger(num_chr_s, num_chr_e, meth)
