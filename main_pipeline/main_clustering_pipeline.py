#!/usr/bin/env python3
# Clusterization of LD matrices with dbscan and hdsbcan
# Updates LOC column in main chr table, prepares .hlist files for hdbscan and dbscan
# Tests haplotypes using plink 1.07
# Consolidate assoc.hap files, sort entries by p-value, wirte csv outfile
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 17, 2023
# Please check README.txt
# python3 main_clustering_pipeline.py -chr 1 -m all -e 0.1 -c 0.1 -ms 2 -mfp /mnt/wd/nsap/in_data1/chr1_matrix.ld -lfp /mnt/wd/nsap/in_data1/chr1.snplist

import pandas as pd
import numpy as np
import datetime
from sklearn.cluster import DBSCAN
from sklearn.metrics.cluster import *
import hdbscan
import sys
import argparse    
import re
from joblib import Memory
import mysql.connector as mysql
from mysql.connector import ProgrammingError
from sqlalchemy import create_engine
import configparser
import os

import SNP_clustering_v092
import hlists_prepare
import assoc_hap_Hdbscan
import merger_assoc_hap_v02

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


def _main_():
  ### Time of start
  time_start = datetime.datetime.now()
  print (f"\nTime of start - {time_start.isoformat(sep=' ', timespec='seconds')}") 
  
  # Parse command line arguments
  parser = argparse.ArgumentParser(description='Parameters for clustering algoritms. Table and file names will also depend on these.')
  parser.add_argument('-chr', '--num_chr', 
                      help = 'number of a chromosome either sole or separated with a \'-\'', 
                      required = True)
  parser.add_argument('-m', '--method', 
                      help = 'possible entries are: hdbscan, dbscan, all', 
                      required = True)
  parser.add_argument('-e', '--eps',
                      help = 'epsilon parameter for DBSCAN',
                      required = True)
  parser.add_argument('-c', '--csm',
                      help = 'cluster_selection_epsilon parameter for HDBSCAN',
                      required = True)
  parser.add_argument('-ms', '--min_samples',
                      help = 'min_samples parameter for DBSCAN and HDBSCAN, min_cluster_size of HDBSCAN also equals to this entry',
                      required = True)
  parser.add_argument('-mfp', '--matrix_file_path',
                      help = 'file path to .ld file',
                      required = True)
  parser.add_argument('-lfp', '--snp_file_path',
                      help = 'file path to .snp_file_path file',
                      required = True)
  args = parser.parse_args()
  
  # Variables
  matrix_file_path = args.matrix_file_path
  snp_file_path = args.snp_file_path
  meth = args.method
  num_chr = args.num_chr
  eps = float(args.eps)
  csm = float(args.csm)
  min_samples = int(args.min_samples)
  
  # Print out parameters
  print('Arguments in use:', '\n', 
      '- matrix filepath', matrix_file_path, '\n',
      '- snplist filepath', snp_file_path, '\n',
      '- method', meth, '\n',
      '- chr', num_chr, '\n',
      '- eps', eps, '\n',
      '- cluster_selection_epsilon', csm, '\n',
      '- min_samples', min_samples, '\n')
  ### Clusterization parameters
  eps_str = str(eps).replace('.', '')
  csm_str = str(csm).replace('.', '')
  ### Defining filename which is also an alias for tables
  filename = f"A_chr{num_chr}_CSM{csm_str}MS{min_samples}_EPS{eps_str}MS{min_samples}"
  # Clustering with SNP_clustering_v0....py
  SNP_clustering_v092.awaiting_input(matrix_file_path,
                                     snp_file_path,
                                     meth,
                                     eps,
                                     csm,
                                     min_samples,
                                     num_chr,
                                     filename)
  # Preparing hlists with hlists_prepare.py
  hlists_prepare.updateTableLOC(num_chr,
                                filename)
  hlists_prepare.hdbscan_hlist(num_chr,
                               filename)
  hlists_prepare.dbscan_hlist(num_chr,
                              filename)
  # Restore and test haplotypes with assoc_hap_Hdbscan.py
  assoc_hap_Hdbscan.assoc_hap(num_chr,
                              filename)
  # Sort and filter to out .csvs with merger_assoc_hap.py
  if meth == 'dbscan' or meth == 'hdbscan':
    merger_assoc_hap_v02.dataframe_merger(num_chr, 
                                          num_chr,                              # the function requires range
                                          meth, 
                                          filename)
  else:
    merger_assoc_hap_v02.dataframe_merger(num_chr, 
                                          num_chr,                              # the function requires range
                                          'dbscan', 
                                          filename)    
    merger_assoc_hap_v02.dataframe_merger(num_chr, 
                                          num_chr,                              # the function requires range
                                          'hdbscan', 
                                          filename) 
  ### Time of end
  time_end = datetime.datetime.now()
  print(f"""\nTime of finish - {time_end.isoformat(sep=' ', timespec='seconds')}. 
            \nTime of executing - {time_end - time_start}.""")
### Create a log file
# with open('main_clustering_pipeline.log', 'w') as lf:
#   print('Logging in the \'main_clustering_pipeline.log\'.')
### Open the log file in append mode
try:
    log_file = open('main_clustering_pipeline.log', 'a')
    ### Create an instance of the custom Logger class
    logger = Logger(log_file)
    ### Set the standard output to the Logger instance
    sys.stdout = logger
    ### _main_()
    _main_()

finally:
  ### Close the log file
  log_file.close()
