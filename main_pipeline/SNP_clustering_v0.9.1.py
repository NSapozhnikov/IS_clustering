#!/usr/bin/env python3
# Clusterization of LD matrices with dbscan and hdsbcan
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: April 06, 2023
# snp_clustering v0.9.1
# Please check README.txt
#        check clusterization parameters and table names
# python3 SNP_clustering_v0.9.1.py /mnt/wd/nsap/in_data1/chr1_matrix.ld /mnt/wd/nsap/in_data1/chr1.snplist all

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

### Clusterization parameters
eps = 0.75
csm = 0.65
eps_str = str(eps).replace('.', '')
csm_str = str(csm).replace('.', '')
min_samples = 2

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

def dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, filename):
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
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_dbscan.clustered", 'w') as fp:
    # print('Writting to the outfile...)
    # for i in range(-1, n_clusters_dbscan):
    #     cluster_indices = np.where(labels == i)[0]
    #     for j in cluster_indices:
    #         fp.write(str(snp_list[j]))
    #         fp.write(' ')
    #         fp.write(str(i))
    #         fp.write('\n')
    ### For output into DB    
    print('\nWritting results to database...')
    ### for splited chromosomes. x+n, where n is number of clusters for other part
    # new_labels = []
    # for x in labels:
    #   if x!=-1:
    #     new_labels.append(x+2048)
    #   else:
    #     new_labels.append(-1)
    # pd_df = pd.DataFrame({'SNP': snp_list,
    #                       'cl_index': new_labels})
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index': labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table1', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE {filename} c
                    JOIN temporary_table1 u ON c.snp = u.snp
                    SET c.dbscan = u.cl_index;"""
    with engine.begin() as conn:
        conn.execute(add_query)
    return n_clusters_dbscan, shil_score_dbscan, CH_dbscan
  
def hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr, filename):
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
    ### For output as file 
    # out_file_path = matrix_file_path.rstrip('.ld')
    # with open(f"{out_file_path}_hdbscan.clustered", 'w') as fp:  
    # print('Writting to the outfile...')
    # for i in range(-1, n_clusters_):
    #     cluster_indeces = np.where(labels == i)[0]
    #     for j in cluster_indeces:
    #         #fp.write(str(snp_list[j]))
    #         #fp.write(' ')
    #         #fp.write(str(i))
    #         #fp.write('\n')
    ### For output into DB 
    print('\nWritting results to database...')
    ### for splited chromosomes. x+n, where n is number of clusters for other part 
    # new_labels = []
    # for x in labels:
    #   if x!=-1:
    #     new_labels.append(x+556)
    #   else:
    #     new_labels.append(-1)
    # pd_df = pd.DataFrame({'SNP': snp_list,
    #                       'cl_index': new_labels})
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
    return n_clusters_hdbscan, shil_score_hdbscan, CH_hdbscan
  
def awaiting_input():
    try:
        matrix_file_path = str(sys.argv[1])
        f = open(matrix_file_path, 'r')
        f.close()
        string_list = matrix_file_path.split('/')
        num_chr = re.findall("\\D*(\\d*)", string_list[-1])[0]
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid matrix file path.')
    try:
        snp_file_path = str(sys.argv[2])
        f = open(snp_file_path, 'r')
        f.close()
        string_list1 = snp_file_path.split('/')
        num_chr1 = re.findall("\\D*(\\d*)", string_list1[-1])[0]
        if num_chr1 != num_chr:
            sys.exit('Chromosome missmatch in file names.')
    except (IndexError, FileNotFoundError, ValueError ):
        sys.exit('Invalid SNP list file path.')    
    try:
        methods = ['dbscan','hdbscan','all']
        typed_input = sys.argv[3].lower()
        print(f"\n{typed_input} algorithm will be used...")
    except IndexError:
        sys.exit(f"Invalid method: {typed_input}.")
    ### Time of start
    time_start = datetime.datetime.now()
    print (f"\nTime of start - {time_start.isoformat(sep=' ', timespec='seconds')}") 
    ### Establishing connection with MySQL
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    engine = create_engine(establish_engine())
    
    ### Defining filename which is also an alias for tables
    filename = f"A_chr{num_chr}_CSM{csm_str}MS{min_samples}_EPS{eps_str}MS{min_samples}"
    
    ### Exctracting files
    try:
        with open (snp_file_path, 'r') as snp_file: 
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
    ################################
    if cursor.fetchone()[0] == 0:
    ################################
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
      n_clusters_dbscan, shil_score_dbscan, CH_dbscan = dbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
    elif typed_input == 'hdbscan':
      print('Reshaping an array...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      corr_matrix = None
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)              
      n_clusters_hdbscan, shil_score_hdbscan, CH_hdbscan = hdbscan_clustering(diss_matrix, snp_list, matrix_file_path, num_chr)
    elif typed_input == 'all':
      print('Reshaping arrays...')
      diss_matrix = 1 - np.abs(corr_matrix, out=corr_matrix)
      corr_matrix = None
      diss_matrix = diss_matrix.reshape(len(snp_list), len(snp_list))
      np.fill_diagonal(diss_matrix, 0)
      n_clusters_dbscan, shil_score_dbscan, CH_dbscan = dbscan_clustering(diss_matrix, 
                                                                          snp_list, 
                                                                          matrix_file_path, 
                                                                          num_chr,
                                                                          filename)
      n_clusters_hdbscan, shil_score_hdbscan, CH_hdbscan = hdbscan_clustering(diss_matrix, 
                                                                              snp_list, 
                                                                              matrix_file_path, 
                                                                              num_chr,
                                                                              filename)
    ### Adding row with scores to the xlsx table 
    # try:
    #   print('\nChecking if scores file exists.')
    #   score_file = open('/mnt/wd/nsap/clustering_scores_CSM041MS5.xlsx', 'r')
    #   score_file.close()
    # except FileNotFoundError:
    #   with open('/mnt/wd/nsap/clustering_scores_CSM041MS5.xlsx', 'w') as fp:
    #     scores_pddf = pd.DataFrame({'CHR':[],
    #                     'Number_of_SNPs':[],
    #                     'DBSCAN_clusters':[],
    #                     'DBSCAN_Shil':[],
    #                     'DBSCAN_CH':[],
    #                     'HDBSCAN_clusters':[],
    #                     'HDBSCAN_Shil':[],
    #                     'HDBSCAN_CH':[]})
    #     scores_pddf.to_excel('/mnt/wd/nsap/clustering_scores_CSM041MS5.xlsx',index=False,engine='openpyxl')
    #     print("File 'clustering_scores_CSM041MS5.xlsx' is created.")
    # scores_file_pddf = pd.read_excel('/mnt/wd/nsap/clustering_scores_CSM041MS5.xlsx',
    #                                   engine='openpyxl')
    # scores_pddf = pd.DataFrame({'CHR':[num_chr],
    #                             'Number_of_SNPs':[len(snp_list)],
    #                             'DBSCAN_clusters':[n_clusters_dbscan],
    #                             'DBSCAN_Shil':[shil_score_dbscan],
    #                             'DBSCAN_CH':[CH_dbscan],
    #                             'HDBSCAN_clusters':[n_clusters_hdbscan],
    #                             'HDBSCAN_Shil':[shil_score_hdbscan],
    #                             'HDBSCAN_CH':[CH_hdbscan]})
    # scores_file_pddf = pd.concat([scores_file_pddf, scores_pddf])
    # scores_file_pddf.to_excel('/mnt/wd/nsap/clustering_scores_CSM041MS5.xlsx',
    #                           index=False,engine='openpyxl')
    # print('Success. Row added.')
    ### Time of end
    time_end = datetime.datetime.now()
    print(f"""\nTime of finish - {time_end.isoformat(sep=' ', timespec='seconds')}. 
              \nTime of executing - {time_end - time_start}.""")

# awaiting_input()   

