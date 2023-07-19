#!/usr/bin/env python3
# Clusterization of LD matrices with dbscan and hdsbcan
# Author: Nikita Sapozhnikov, info@inzilico.com
# Date: May 29, 2023
# best_params_search v0.2
# Please check README.txt
# python3 best_params_search.py
# for some reason HDBSCAN doesn't accept np.arange() list
# used parameters are:
# params = {'eps': [0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 ,
#                   0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95],
#           'min_samples': [ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 
#                            16, 17, 18, 19]}
# total count of rows in the consolidated file: amount of chromosomes*len(eps)*len(min_samples)*2 clustering methods

import numpy as np
import datetime
from sklearn.cluster import  DBSCAN 
from sklearn.metrics.cluster import *
import hdbscan
import pandas as pd
import sys
import argparse    
from getpass import getpass
import re
import configparser
from joblib import Memory
from sklearn.model_selection import GridSearchCV

### Define a custom class to duplicate output to both console and log file
class Logger(object):
    def __init__(self, log_file):
        self.terminal = sys.stdout
        self.log_file = log_file

    def write(self, message):
        self.terminal.write(message)
        self.log_file.write(message)
    
    def flush(self):
        self.terminal.flush()
        if self.log_file and not self.log_file.closed:
            self.log_file.flush()

    def close(self):
        self.flush()
        if self.log_file and not self.log_file.closed:
            self.log_file.close()

try:
  log_file = open('best_params_clustering.log', 'a')
  ### Create an instance of the custom Logger class
  logger = Logger(log_file)
  ### Set the standard output to the Logger instance
  sys.stdout = logger

  main_fp = '/mnt/wd/nsap/Clustering/Parameters_search/'
  
  ### Initialize an empty DataFrame for the results
  all_chr_df = pd.DataFrame(columns=['chr','algorithm','clusters','eps','min_samples','sil_score','CH_score'])
  
  ### Define the parameter grid
  params = {'eps': [0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 ,
                    0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95],
            'min_samples': [ 2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]}
  
  ### check if the out-file .xlsx exists
  try:
    with open(main_fp+'best_params_clustering_A.xlsx', 'r') as xlsx_file:
      print('"best_params_clustering_A.xlsx" already exits. It will be overwrited!')
  except FileNotFoundError:
    with open(main_fp+'best_params_clustering_A.xlsx', 'w') as xlsx_file:
      print('"best_params_clustering_A.xlsx" doesn\'t exit. It has been created.')
          
  
  for num_chr in range(1,23):
      ### Initialize an empty DataFrame for the results
      data_df = pd.DataFrame(columns=['chr','algorithm','clusters','eps','min_samples','sil_score','CH_score'])
      ### open in-files with snplist and matrix for a chromosome
      snp_file_path = f"/mnt/wd/nsap/in_data1/chr{num_chr}.snplist"
      matrix_file_path = f"/mnt/wd/nsap/in_data1/chr{num_chr}_matrix.ld"
      try:
          with open(snp_file_path, 'r') as snp_file:
              print('\nExtracting SNP list data from file...')
              snp_list = np.loadtxt(snp_file, dtype=str)
              print(f"The SNP list for {str(num_chr)} chromosome is: \n{snp_list}")
              print(f"With length of: {len(snp_list)}")
      except FileNotFoundError:
          sys.exit('No such file.')      
      try:
          with open(matrix_file_path, 'r') as corr_file:
              print('Extracting matrix data from file...')
              corr_matrix = np.fromfile(corr_file, sep=' ')
      except FileNotFoundError:
          sys.exit('No such file.')
      ### transforming into a dissimilarity matrix    
      np.nan_to_num(corr_matrix, copy=False)  
      print('Reshaping an array...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)
      
      ### Initialize mem instance for hdbscan
      mem = Memory(location='./cachedir')
      ### Iterate over the parameter grid
      for eps in params['eps']:
          for min_samples in params['min_samples']:
              print('eps: ', eps)
              print('min_samples: ', min_samples)
              ### DBSCAN estimator
              dbscan = DBSCAN(eps=eps, 
                              min_samples=min_samples, 
                              metric='precomputed', 
                              n_jobs=10).fit(diss_matrix) 
              ### Calculate scores that measure the quality of the clustering
              labels = dbscan.labels_
              n_clusters_dbscan = len(set(labels)) - (1 if -1 in labels else 0)
              print('Estimated number of DBSCAN clusters: %d' % n_clusters_dbscan)
              dbscan_sil_score = silhouette_score(diss_matrix, labels, metric='precomputed') if len(set(labels)) > 1 else 0
              dbscan_CH_score = calinski_harabasz_score(diss_matrix, labels) if len(set(labels)) > 1 else 0
              ### append scores to a DataFrame
              data_df.loc[len(data_df.index)] = [num_chr,
                                                 'dbscan',
                                                 n_clusters_dbscan,
                                                 eps,
                                                 min_samples,
                                                 dbscan_sil_score,
                                                 dbscan_CH_score]
              print(data_df.loc[len(data_df.index) - 1])
              ### HDSBCAN estimator
              clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples, 
                                          min_samples=min_samples,
                                          cluster_selection_epsilon=eps,
                                          metric='precomputed', 
                                          cluster_selection_method='leaf',
                                          core_dist_n_jobs=10, 
                                          memory=mem).fit(diss_matrix)
              ### Calculate scores that measure the quality of the clustering
              labels = clusterer.labels_
              n_clusters_hdbscan = len(set(labels)) - (1 if -1 in labels else 0)
              print('Estimated number of HDBSCAN clusters: %d' % n_clusters_hdbscan)
              hdbscan_sil_score = silhouette_score(diss_matrix, labels, metric='precomputed') if len(set(labels)) > 1 else 0
              hdbscan_CH_score = calinski_harabasz_score(diss_matrix, labels) if len(set(labels)) > 1 else 0
              ### append scores to list
              data_df.loc[len(data_df.index)] = [num_chr,
                                                 'hdbscan',
                                                 n_clusters_hdbscan,
                                                 eps,
                                                 min_samples,
                                                 hdbscan_sil_score,
                                                 hdbscan_CH_score]
              print(data_df.loc[len(data_df.index) - 1])
      ### Write an out-file for the chromosome
      print(f"Writing an out-file for chromosome {num_chr}...")
      data_df.to_excel(main_fp+f"best_params_clustering_A_chr{num_chr}.xlsx",
                       sheet_name='Sheet1',
                       index=False,
                       header=True)
      all_chr_df = pd.concat([all_chr_df, data_df], ignore_index=True)
      ### clear mem instance
      mem.clear(warn=False)
  print(all_chr_df)
  ### Write a consolidated out-file
  print('Writing a consolidated out-file...')
  all_chr_df.to_excel(main_fp+'best_params_clustering_A.xlsx', 
                   sheet_name='Sheet1', 
                   index=False, 
                   header=True)
  print('Done!')
except Exception:
    traceback.print_exc()
finally:
  ### Close the log file
  log_file.close()
        
        
        
        
