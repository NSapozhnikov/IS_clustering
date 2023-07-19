#!/usr/bin/env python3
# Automation of main_clustering_pipeline.py
# Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
# Date: July 17, 2023
# Please check README.txt
# python3 all_chr_pipe.py

import os

for i in range(1,23):
  cmd = f"python3 main_clustering_pipeline.py -chr {i} -m all -e 0.75 -c 0.7 -ms 2 -mfp /mnt/wd/nsap/in_data1/chr1_matrix.ld -lfp /mnt/wd/nsap/in_data1/chr1.snp_file_path"
  print(cmd)
  os.system(cmd)
