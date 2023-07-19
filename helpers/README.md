### This folder contains helper scripts

# make_bed_imp1.py
  - Make plink binary format files from merged using plink
  python3 make_bed_imp1.py

# SNP_LOC.py
  - Creates/Updates tables with SNPs locations in MySQL db
  python3 SNP_LOC.py
  
# CHR_START_check.py
  - Reads .csv files after merger_assoc_hap.py and prints out number of clusters
  which are unique CHR:START combinations. 
  - This script is for comparing number of clusters after clusterization and 
  after forming and testing haplotypes as Plink does not test haplotypes with 
  length 1, however there are clusters with only 1 SNP

# count_clusters_from_DB.py
  - Reads MySQL tables and prints out number of clusters right after clusterization
