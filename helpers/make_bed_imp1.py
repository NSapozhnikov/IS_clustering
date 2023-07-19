#!/usr/bin/env python3
# Make plink binary format files from merged using plink 
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 11, 2023
# make_bed_imp1
# Please check README.txt
# python3 /mnt/wd/nsap/Clustering/helpers/make_bed_imp1.py

import os

fp = '/mnt/wd/nsap/'

for num_chr in range(1,23):
    os.system(f"/mnt/wd/nsap/./plink --bfile {fp}imp1/merged --chr {num_chr} --make-bed --out {fp}imp1/chr{num_chr}")
    print(f"Element {num_chr} has been plinked...")
    os.remove(f"{fp}imp1/chr{num_chr}.log")
    try:
      os.remove(f"{fp}imp1/chr{num_chr}.nosex")
    except FileNotFoundError:
      print('No .nosex file found.')
print('Done!')
