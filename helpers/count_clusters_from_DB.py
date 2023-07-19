import numpy as np
import pandas as pd
import mysql.connector as mysql
import configparser

def connectToBD_MySQL():
    config = configparser.ConfigParser()
    config.read('credentials.ini')  
    return mysql.connect(
            host=config['DEFAULT']['host'],
            user=config['DEFAULT']['user'],
            password=config['DEFAULT']['password'],
            database='cl')
            
counter = 0

### Establishing connection with MySQL
connection = connectToBD_MySQL()
cursor = connection.cursor()
for num_chr in range(1,23):
  if num_chr == 6:
    read_query = f"""SELECT dbscan FROM clusters_chr6_CSM041MS5_noMHC
                     ORDER BY
                      CAST(dbscan as decimal) DESC
                     LIMIT 1;"""
  else:
    read_query = f"""SELECT dbscan FROM clusters_chr{num_chr}_CSM041MS5
                     ORDER BY
                      CAST(dbscan as decimal) DESC
                     LIMIT 1;"""
                     
  cursor.execute(read_query)
  # cluster numbers count 0-origin
  clusters = int(cursor.fetchall()[0][0]) + 1 
  connection.commit()
  print(f"Chromosome {num_chr} DBSCAN clustering has {clusters} clusters.")
  counter += clusters
print(f"{counter} clusters in totall.")

