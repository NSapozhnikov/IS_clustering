#!/usr/bin/env python3
# Clusterization of LD matrices with dbscan and hdsbcan
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 17, 2023
# snp_clustering v0.9.2 - pipeline edition
# Please check README.txt
#        check clusterization parameters and table names

import numpy as np
import datetime
from sklearn.cluster import DBSCAN
from sklearn.metrics.cluster import *
import hdbscan
import pandas as pd
import sys
import argparse    
import mysql.connector as mysql
import re
import configparser
from sqlalchemy import create_engine
from joblib import Memory

def connectToBD_MySQL():
    config = configparser.ConfigParser()
    config.read('credentials.ini')  
    return mysql.connect(
            host=config['DEFAULT']['host'],
            user=config['DEFAULT']['user'],
            password=config['DEFAULT']['password'],
            database='cl')

def establish_engine():
    config = configparser.ConfigParser()
    config.read('credentials.ini')   
    engine = config['DEFAULT']['engine']
    return engine

# DBSCAN
def dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, filename, eps, min_samples):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with DBSCAN...')
    db = DBSCAN(eps=eps, 
                min_samples=min_samples, 
                metric='precomputed', 
                n_jobs=6).fit(diss_matrix)
    labels = db.labels_
    n_clusters_dbscan = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_dbscan)
    ### Scoring
    try:
        shil_score_dbscan = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', shil_score_dbscan)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH_dbscan = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', CH_dbscan)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output into DB    
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index': labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table1', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE {filename} c
                    JOIN temporary_table1 u ON c.snp = u.snp
                    SET c.dbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)

# HDBSCAN
def hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, filename, csm, min_samples):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with HDBSCAN...')
    mem = Memory(location='./cachedir')
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples, 
                                min_samples=min_samples,
                                cluster_selection_epsilon=csm,
                                metric='precomputed', 
                                cluster_selection_method='leaf',
                                core_dist_n_jobs=6, memory=mem).fit(diss_matrix)
    labels = clusterer.labels_
    n_clusters_hdbscan = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_hdbscan)
    ### Scoring
    try:
        shil_score_hdbscan = silhouette_score(diss_matrix, labels, metric='precomputed') 
        print('Silhouette score: ', shil_score_hdbscan)
    except:
        print('Failed to compute Silhouette score.')
    try:
        CH_hdbscan = calinski_harabasz_score(diss_matrix, labels)
        print('Calinski-Harabasz score: ', CH_hdbscan)
    except:
        print('Failed to compute Calinski-Harabasz score.')
    ### For output into DB 
    print('\nWritting results to database...')
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index': labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table1', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE {filename} c
                    JOIN temporary_table1 u ON c.snp = u.snp
                    SET c.hdbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)
    mem.clear(warn=False)

def awaiting_input(matrix_file_path, snp_file_path, meth, eps, csm, min_samples, num_chr, filename):
    # Check if matrix filepath is correct
    try:
        with open(matrix_file_path, 'r') as fp:
          string_list = matrix_file_path.split('/')
          num_chr_ = re.findall("\\D*(\\d*)", string_list[-1])[0]
          # Check num_chr
          if int(num_chr) != int(num_chr_):
            print(f"Warning! Entered parameter \'{num_chr}\' not equals to infered from the file path \'{num_chr_}\'!")
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid matrix file path.')
    # Check if SNPlist filepath is correct
    try:
        with open(snp_file_path, 'r') as fp:
          string_list1 = snp_file_path.split('/')
          num_chr1 = re.findall("\\D*(\\d*)", string_list1[-1])[0]
          if num_chr1 != num_chr:
              sys.exit('Chromosome missmatch in file names.')
    except (IndexError, FileNotFoundError, ValueError) as error:
        # sys.exit('Invalid SNP list file path.')    
        print(error)
    try:
        methods = ['dbscan','hdbscan','all']
        typed_input = meth.lower()
        print(f"\n{typed_input} algorithm will be used...")
    except IndexError:
        sys.exit(f"Invalid method: {typed_input}.")

    ### Establishing connection with MySQL
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    engine = create_engine(establish_engine())

    ### Exctracting files
    with open (snp_file_path, 'r') as snp_file: 
        print('\nExtracting SNP list data from file...')
        snp_list = np.loadtxt(snp_file, dtype=str)        
        print(f"The SNP list for {str(num_chr)} chromosome is: \n{snp_list}")
        print(f"With length of: {len(snp_list)}")
    with open(matrix_file_path, 'r') as corr_file:
        print('Extracting matrix data from file...')
        corr_matrix = np.fromfile(corr_file, sep=' ')
    np.nan_to_num(corr_matrix, copy=False)
    diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
    diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
    np.fill_diagonal(diss_matrix, 0)
    corr_matrix = corr_matrix.reshape(len(snp_list), len(snp_list))
    np.fill_diagonal(corr_matrix, 1)
    
    ### Check if table exists in the schema, create if not  
    cursor.execute(f"""SELECT * 
                       FROM information_schema.tables
                       WHERE table_schema = 'cl' 
                       AND table_name = '{filename}'
                       LIMIT 1;""")
    if cursor.fetchone() == None:
        print(f"No table in the schema was found for chromosome {num_chr}.", 
               'A new table will be created...', sep='\n')
        cursor.execute(f"""CREATE TABLE {filename}(
                                SNP VARCHAR(32)
                                ,LOC VARCHAR(32)
                                ,dbscan VARCHAR(32)
                                ,hdbscan VARCHAR(32));""")
        connection.commit()
    ### Fill SNP column with rs-IDs from snplist file
    cursor.execute(f"""SELECT count(*) AS total FROM {filename};""")
    if cursor.fetchone()[0] == 0:
      baseSQL = f"INSERT INTO cl.{filename} (SNP) VALUES"
      for snp in snp_list:
          baseSQL += "('" + snp + "'),"
      baseSQL = baseSQL.rstrip(',')
      with engine.begin() as conn:
          conn.execute(baseSQL)

    if typed_input == 'dbscan':
      print('Reshaping an array...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      corr_matrix = None
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)
      dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, eps, min_samples)
    elif typed_input == 'hdbscan':
      print('Reshaping an array...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      corr_matrix = None
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)              
      hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, csm, min_samples)
    elif typed_input == 'all':
      print('Reshaping arrays...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      corr_matrix = None
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)
      dbscan_clustering(diss_matrix,snp_list,matrix_file_path,num_chr,filename, eps, min_samples)
      hdbscan_clustering(diss_matrix,snp_list,matrix_file_path,num_chr,filename, csm, min_samples)
