#!/usr/bin/env python3
# Clusterization of LD matrices with dbscan and hdsbcan
# Author: Nikita Sapozhnikov, info@inzilico.com
# Date: June 7, 2023
# best_params_clustering v0.1
# Please check README.txt
# python3 best_params_clustering_1.py

import numpy as np
import datetime
from sklearn.cluster import  DBSCAN 
from sklearn.metrics.cluster import *
import hdbscan
import pandas as pd
import sys
import configparser
from joblib import Memory
import seaborn as sns
import matplotlib.pyplot as plt
from sqlalchemy import create_engine
import mysql.connector as mysql

main_fp = '/mnt/wd/nsap/'

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

### MySQL connector
def connectToBD_MySQL():
    config = configparser.ConfigParser()
    config.read('credentials.ini')  
    return mysql.connect(
            host=config['DEFAULT']['host'],
            user=config['DEFAULT']['user'],
            password=config['DEFAULT']['password'],
            database='cl')
            
### MySQL dump engine
def establish_engine():
    config = configparser.ConfigParser()
    config.read('credentials.ini')   
    engine = config['DEFAULT']['engine']
    return engine

### DBSCAN clusterer
def dbscan_clustering(diss_matrix, snp_list, num_chr, eps, min_samples):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with DBSCAN...')
    db = DBSCAN(eps=eps, 
                min_samples=min_samples, 
                metric='precomputed', 
                n_jobs=10).fit(diss_matrix)
    labels = db.labels_
    n_clusters_dbscan = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_dbscan)
    ### Scoring
    try:
      shil_score_dbscan = silhouette_score(diss_matrix, labels, metric='precomputed')
      print('Silhouette score: ', shil_score_dbscan)
    except:
        print('Failed to compute Silhouette score.')
        shil_score_dbscan = 0
    try:
      CH_dbscan = calinski_harabasz_score(diss_matrix, labels)
      print('Calinski-Harabasz score: ', CH_dbscan)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        calinski_harabasz_score = 0
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index': labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table1', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE imp1_chr{num_chr}_best_sil c
                    JOIN temporary_table1 u ON c.snp = u.snp
                    SET c.dbscan = u.cl_index;"""
    print(f"Updating MySQL table for {num_chr} chromosome...")
    with engine.begin() as conn:
        conn.execute(add_query)
    print('The table has been updated.')
    return n_clusters_dbscan, shil_score_dbscan, CH_dbscan

### HDBSCAN clusterer
def hdbscan_clustering(diss_matrix, snp_list, num_chr, eps, min_samples):
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    print('\nPerforming clustering with HDBSCAN...')
    mem = Memory(location='./cachedir')
    clusterer = hdbscan.HDBSCAN(min_cluster_size=8, 
                                min_samples=min_samples,
                                cluster_selection_epsilon=eps,
                                metric='precomputed', 
                                cluster_selection_method='leaf',
                                core_dist_n_jobs=10, 
                                memory=mem).fit(diss_matrix)
    labels = clusterer.labels_
    n_clusters_hdbscan = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_hdbscan)
    ### Scoring
    try:
      shil_score_hdbscan = silhouette_score(diss_matrix, labels, metric='precomputed')
      print('Silhouette score: ', shil_score_hdbscan)
    except:
        print('Failed to compute Silhouette score.')
        shil_score_hdbscan = 0
    try:
      CH_hdbscan = calinski_harabasz_score(diss_matrix, labels)
      print('Calinski-Harabasz score: ', CH_hdbscan)
    except:
        print('Failed to compute Calinski-Harabasz score.')
        CH_hdbscan = 0
    pd_df = pd.DataFrame({'SNP': snp_list,
                          'cl_index': labels})
    engine = create_engine(establish_engine())
    pd_df.to_sql(name='temporary_table1', con=engine, schema='cl', if_exists='replace', index=False)
    add_query = f"""UPDATE imp1_chr{num_chr}_best_sil c
                    JOIN temporary_table1 u ON c.snp = u.snp
                    SET c.hdbscan = u.cl_index;"""
    print(f"Updating MySQL table for {num_chr} chromosome...")
    with engine.begin() as conn:
        conn.execute(add_query)
    mem.clear(warn=False)
    print('The table has been updated.')
    return n_clusters_hdbscan, shil_score_hdbscan, CH_hdbscan

### Plot heatmaps
def plot_HM(meth, num_chr, sub_df):
    ### Plot the heatmap for Silhouette Score
    ### Create a figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    ### Pivot a heatmap
    sil_score_HM = sub_df.pivot_table(index='eps', 
                                      columns='min_samples', 
                                      values='sil_score')
    clusters_HM = sub_df.pivot_table(index='eps', 
                                     columns='min_samples', 
                                     values='clusters')
    ### construct a mask for (hide) values where clusters == 0
    mask_shape = sil_score_HM.shape
    mask = np.array(sub_df['clusters'] == 0).reshape(mask_shape)
    ### Plot the Silhouette Score heatmap
    sns.heatmap(sil_score_HM, 
                cmap='viridis', 
                annot=clusters_HM, 
                linewidths=0.1,
                fmt=".0f",
                mask=mask,
                ax=axes[0])
    axes[0].set_title('Silhouette Score')
    axes[0].set_xlabel('min_samples')
    axes[0].set_ylabel('eps')    
    ### Pivot a heatmap
    CH_score_HM = sub_df.pivot_table(index='eps', 
                                     columns='min_samples', 
                                     values='CH_score')
    ### Plot the CH Score heatmap
    sns.heatmap(CH_score_HM, 
                cmap='viridis', 
                annot=clusters_HM, 
                linewidths=0.1,
                fmt=".0f", 
                mask=mask,
                ax=axes[1])
    axes[1].set_title('CH Score')
    axes[1].set_xlabel('min_samples')
    axes[1].set_ylabel('eps')
    plt.tight_layout()
    # Add a title to the composition
    fig.suptitle(f"Chromosome {num_chr} {meth.upper()} scoring", fontsize=16, fontweight='bold',y=0.99)
    ### Save the figure as a PDF file
    plt.savefig(f"{main_fp}Clustering/score_heatmaps/{num_chr}_{meth}_HM.pdf", format='pdf')
    plt.show()

### Prepare in-matrix and MySQL tables
def before_clustering(num_chr):
    ### Clustering with best params for sil_score
    connection = connectToBD_MySQL()
    cursor = connection.cursor()
    engine = create_engine(establish_engine())
    ### Exctracting files
    snp_file_path = f"{main_fp}in_data1/chr{num_chr}.snplist"
    matrix_file_path = f"{main_fp}in_data1/chr{num_chr}_matrix.ld"
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
    ### Check if table exists in the schema, create if not
    cursor.execute(f"""SELECT *
                       FROM information_schema.tables
                       WHERE table_schema = 'cl'
                       AND table_name = 'imp1_chr{num_chr}_best_sil'
                       LIMIT 1;""")
    if cursor.fetchone() == None:
        print(f"No table in the schema was found for chromosome {num_chr}.",
               'A new table will be created...', sep='\n')
        cursor.execute(f"""CREATE TABLE imp1_chr{num_chr}_best_sil(
                                SNP VARCHAR(32)
                                ,LOC VARCHAR(32)
                                ,dbscan VARCHAR(32)
                                ,hdbscan VARCHAR(32));""")
        connection.commit()
    ### Fill SNP column with rs-IDs from snplist file
    cursor.execute(f"""SELECT count(*) AS total FROM imp1_chr{num_chr}_best_sil;""")
    ################################ comment this if appending SNPs of the second half of the splitted chromosome
    if cursor.fetchone()[0] == 0:
    ################################
      baseSQL = f"INSERT INTO cl.imp1_chr{num_chr}_best_sil (SNP) VALUES"
      for snp in snp_list:
          baseSQL += "('" + snp + "'),"
      baseSQL = baseSQL.rstrip(',')
      with engine.begin() as conn:
          conn.execute(baseSQL)
    return diss_matrix, snp_list

### function to compare score from an in-file and score from clustering
def compare_scores(a,b):
    if a == b:
      print('Scores are equal.')
    elif a > b:
      print(f"Chosen from an in-file score ({a}) is more than score from the clustering ({b}).")
    elif b > a:
      print(f"Chosen from an in-file score ({a}) is less than score from the clustering ({b}).")

### Create a log file
with open('best_params_clustering_1.log', 'w') as lf:
  print('Logging in the \'best_params_clustering_1.log\'.')
### Open the log file in append mode
try:
    log_file = open('best_params_clustering_1.log', 'a')
    ### Create an instance of the custom Logger class
    logger = Logger(log_file)
    ### Set the standard output to the Logger instance
    sys.stdout = logger
  
    ### read scores file
    print('Reading consolidated scores file...')
    pd_df = pd.read_excel(main_fp+'helpers/best_params_clustering_A.xlsx')
    
    ### for each chromosome 
    for num_chr in range(1,23):
        print(f"\n\nWorking with chromosome {num_chr}...")
        ### Prepare in-matrix and MySQL tables
        diss_matrix, snp_list = before_clustering(num_chr)
        ### for each method find max scores
        for meth in set(pd_df['algorithm']):
            print(f"\n{meth.upper()}")
            cond1 = pd_df['chr'] == num_chr
            cond2 = pd_df['algorithm'] == meth
            sub_df = pd_df.loc[cond1 & cond2]
            sil_score_max = sub_df['sil_score'].max()
            print('The best Silhouette Score is: ', sil_score_max)
            eps_sil = sub_df.loc[sub_df['sil_score'] == sil_score_max, 'eps']
            CH_score_max = sub_df['CH_score'].max()
            print('The best Calinski-Harabasz score is: ', CH_score_max)
            eps_CH = sub_df.loc[sub_df['CH_score'] == CH_score_max, 'eps']
            ### marking the best parameters for clustering
            ### if several values for the best param, take the middle one
            ### e.g. [0,1,2,3] -> 2
            eps = list(set(eps_sil))[int(len(eps_sil) / 2)]
            min_samples = list(set(sub_df.loc[sub_df['sil_score'] == sil_score_max, 'min_samples']))
            min_samples = min_samples[int(len(min_samples) / 2)]
            print(f"Best chosen parameters are: \teps: {eps} \tmin_samples: {min_samples}")
            
            
            ### Plot heatmaps
            # plot_HM(meth, num_chr, sub_df)
            
            
            ### Clusterize with  best params
            if meth == 'dbscan':
              ### DBSCAN clustering
              n_clusters_dbscan, shil_score_dbscan, CH_dbscan = dbscan_clustering(diss_matrix, snp_list, num_chr, eps, min_samples)
              ### Check if scores from in-file table are the same for selected parameters
              compare_scores(sil_score_max, shil_score_dbscan)
              compare_scores(CH_score_max, CH_dbscan)
            elif meth == 'hdbscan':
              ### HDBSCAN clustering
              n_clusters_hdbscan, shil_score_hdbscan, CH_hdbscan = hdbscan_clustering(diss_matrix, snp_list, num_chr, eps, min_samples)
              ### Check if scores from in-file table are the same for selected parameters
              compare_scores(sil_score_max, shil_score_hdbscan)
              compare_scores(CH_score_max, CH_hdbscan)
        ### Adding row with scores to the xlsx table
        try:
            ### check if scores file
            print('\nChecking if scores file exists.')
            score_file = open(main_fp + 'clustering_scores_imp1.xlsx', 'r')
            score_file.close()
        except FileNotFoundError:
            ### create scores file if there is none
            with open(main_fp + 'clustering_scores_imp1.xlsx', 'w') as fp:
                scores_pddf = pd.DataFrame({'CHR':[],
                                'Number_of_SNPs':[],
                                'DBSCAN_clusters':[],
                                'DBSCAN_Shil':[],
                                'DBSCAN_CH':[],
                                'HDBSCAN_clusters':[],
                                'HDBSCAN_Shil':[],
                                'HDBSCAN_CH':[]})
                scores_pddf.to_excel(main_fp + 'clustering_scores_imp1.xlsx',
                                     index=False,
                                     engine='openpyxl')
                print("File \'clustering_scores_imp1.xlsx\' is created.")
        scores_file_pddf = pd.read_excel(main_fp + 'clustering_scores_imp1.xlsx',
                                          engine='openpyxl')
        scores_pddf = pd.DataFrame({'CHR':[num_chr],
                                    'Number_of_SNPs':[len(snp_list)],
                                    'DBSCAN_clusters':[n_clusters_dbscan],
                                    'DBSCAN_Shil':[shil_score_dbscan],
                                    'DBSCAN_CH':[CH_dbscan],
                                    'HDBSCAN_clusters':[n_clusters_hdbscan],
                                    'HDBSCAN_Shil':[shil_score_hdbscan],
                                    'HDBSCAN_CH':[CH_hdbscan]})
        scores_file_pddf = pd.concat([scores_file_pddf, scores_pddf])
        scores_file_pddf.to_excel(main_fp + 'clustering_scores_imp1.xlsx',
                                  index=False,
                                  engine='openpyxl')
        print('Scores row has been added.')
        break

finally:
  ### Close the log file
  log_file.close()
