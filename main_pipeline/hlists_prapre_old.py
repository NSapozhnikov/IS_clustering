#!/usr/bin/env python3
# Updates LOC column in main chr table, prepares .hlist files for hdbscan and dbscan
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 11, 2023
# hlists_prepare
# Please check README.txt
#        check file and table names
# python3 hlists_prepare.py

import pandas as pd
import mysql.connector as mysql
from mysql.connector import ProgrammingError
from sqlalchemy import create_engine
import configparser

main_fp = '/mnt/wd/nsap/Clustering/main_pipeline/'

### Clusterization parameters
# eps = 0.75
# eps_str = str(eps).replace('.', '')
# min_samples = 2

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
  
connection = connectToBD_MySQL()
cursor = connection.cursor()
engine = create_engine(establish_engine())

def updateTableLOC(num_chr):
    ### Update LOC column
    try:
      upd_query = f"""UPDATE {filename} v
                      JOIN SNP_LOCATIONS_imp1_{num_chr} t ON v.SNP = t.SNP
                      SET v.LOC = t.LOC;"""
      with engine.begin() as conn:
          print(f"Joinning SNP locations for '{filename}' table...")
          conn.execute(upd_query)
    except ProgrammingError:
      print(f"There is no table for chromosome {num_chr} in the schema. Continuing...")
    
def hdbscan_hlist(num_chr):
    ### For HDBSCAN
    try:
      read_query = f"""SELECT SNP, hdbscan, LOC
                       FROM {filename}
                       WHERE (hdbscan NOT LIKE -1)
                       ORDER BY CAST(hdbscan AS DECIMAL) ASC,
                                CAST(LOC AS DECIMAL) ASC;"""
      cursor.execute(read_query)
      clusters__list = cursor.fetchall()
      connection.commit()
      ### Create a dict with unique cluster labels
      counter = {}
      for i in range(len(clusters__list)):
        if clusters__list[i][2] != None:
          key = clusters__list[i][1]
          ## Append (SNP, LOC) tuple by unique keys
          if not key in counter:
              counter[key] = []
          counter[key].append([clusters__list[i][0], clusters__list[i][2]])
      ### Write to the outfile all values with unique keys with more than 1 element
      ### with hlist format for plink's 1.07 wildcard specification --hap flag
      with open(f"{main_fp}hlists/{filename}_hdbscan.hlist", 'w') as fp:
          clusters = 0
          for key in counter:
              if len(counter[key]) > 1:
                  meta_list = counter[key]
                  clusters+=1
                  for j in range(len(meta_list)-1):
                      fp.write(f"** {num_chr}:")
                      snp = f"{meta_list[0][1]}-{meta_list[-1][1]}"
                      fp.write(snp)
                      for v in range(len(meta_list)):
                          fp.write(' ')
                          fp.write(meta_list[v][0])
                      fp.write('\n')
                      meta_list.pop()
          print(f"{clusters} HDBSCAN clusters for chromosome {num_chr} has been written to the outfile.")
    except ProgrammingError:
      print(f"There is no table for chromosome {num_chr} in the schema. Continuing...")

def dbscan_hlist(num_chr):    
    ### For DBSCAN    
    try:
      read_query = f"""SELECT SNP, dbscan, LOC
                       FROM {filename}
                       WHERE (dbscan NOT LIKE -1)
                       ORDER BY CAST(dbscan AS DECIMAL) ASC,
                                CAST(LOC AS DECIMAL) ASC;"""
      cursor.execute(read_query)
      clusters__list = cursor.fetchall()
      connection.commit()
      ### Create a dict with unique cluster labels
      counter = {}
      for i in range(len(clusters__list)):
        if clusters__list[i][2] != None:
          key = clusters__list[i][1]
          ## Append (SNP, LOC) tuple by unique keys
          if not key in counter:
              counter[key] = []
          counter[key].append([clusters__list[i][0], clusters__list[i][2]])
      ### Write to the outfile all values with unique keys with more than 1 element
      ### with hlist format for plink's 1.07 wildcard specification --hap flag
      with open(f"{main_fp}hlists/{filename}_dbscan.hlist", 'w') as fp:
          clusters = 0
          for key in counter:
              if len(counter[key]) > 1:
                  meta_list = counter[key]
                  clusters+=1
                  for j in range(len(meta_list)-1):
                      fp.write(f"** {num_chr}:")
                      snp = f"{meta_list[0][1]}-{meta_list[-1][1]}"
                      fp.write(snp)
                      for v in range(len(meta_list)):
                          fp.write(' ')
                          fp.write(meta_list[v][0])
                      fp.write('\n')
                      meta_list.pop()
          print(f"{clusters} DBSCAN clusters for chromosome {num_chr} has been written to the outfile.")
    except ProgrammingError:
      print(f"There is no table for chromosome {num_chr} in the schema. Continuing...")

# for num_chr in range(1,23):
num_chr = input('num_chr: ')
filename = f"clusters_chr1_sq"
updateTableLOC(num_chr)
hdbscan_hlist(num_chr)
dbscan_hlist(num_chr)
