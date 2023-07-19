#!/usr/bin/env python3
# Creates/Updates tables with SNPs locations in MySQL db
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 11, 2023
# SNP_LOC
# Please check README.txt
# python3 SNP_LOC.py

import pandas as pd
import mysql.connector as mysql
from sqlalchemy import create_engine

def establish_engine():
    config = configparser.ConfigParser()
    config.read('credentials.ini')   
    engine = config['DEFAULT']['engine']
    return engine
  
engine = create_engine(establish_engine())
for num_chr in range(1,23):
  with open(f"/mnt/wd/nsap/imp1/chr{num_chr}.bim", 'r') as fp:
    pd_df = pd.read_csv(fp, header=None, sep='\t', lineterminator='\n', 
                        names=['CHR', 'SNP', 'a', 'LOC', 'b', 'c'])
    pd_df = pd_df.drop(columns=['a','b','c'])
    print(pd_df)
    pd_df.to_sql(name=f"SNP_LOCATIONS_imp1_{num_chr}", con=engine, schema='cl', 
                 if_exists='replace', index=False)
    
