# IS_clustering
A pipeline for clustering squared LD matricies


### Main pipeline: clustering -> hlists preparation -> haplotype testing -> consolidationg haplotypes after Bonferroni


# 1) Squared matrices preparation (-> /mnt/wd/nsap/in_data1 or in_data2)
     python3 plink_matrix_snplist_imp1.py
     

# 2) Chromosome clusterization (SNP_clustering_v0.9.1.py)
     - SNP_clustering_v092.py is a newer version streamlined for the pipeline. v091
       is for the singular usage.
     python3 SNP_clustering_v091.py \
             /mnt/wd/nsap/in_data1/chr1_matrix.ld \
             /mnt/wd/nsap/in_data1/chr1.snplist all


# 3) Writing hlist files (hlists_prepare.py -> hlists/)
     python3 hlists_prepare.py  
   
   
# 4) Haplotype testing (/mnt/wd/nsap/plink-1.07-x86_64/assoc_hap_Hdbscan.py -> hap_assoc/)
     - This step uses plink 1.07 for haplotype testing, thus this 
       automatization script is located in plink's folder
     python3 assoc_hap_Hdbscan.py


# 5) Consolidate and sort by p-value (merger_assoc_hap.py -> hap_assoc/csvs/)
     - merger_assoc_hap_v01.py is an older version. v02 can be used safely.
     python3 merger_assoc_hap.py -n 1 -m all -f A_chr1_CSM07MS2_EPS75MS2
     
     
# 6a) Sum up amount of clusters
     - pipeline_csvs_merger.py is used to merge csvs for different chromosomes. It uses a
       predefined prefix. Correct the prefix for different data sets.
     python3 pipeline_csvs_merger.py


# 6b) Sum up amount of clusters and perform bonferroni correction
     - proc-omnibus-03.R is used to merge csvs for different chromosomes. It uses a
       predefined prefix. Correct the prefix for different data sets. It also performs
       Bonferroni correction and writes an outfile. Plots p-value counts histogram.
     - requires helpers.R

# To automatically perform steps 1-5 on all chromosomes use all_chr_pipe.py. But keep in mind the clusterization parameters.
     python3 all_chr_pipe.py


### Helper scripts

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


### Parameters search (scripts for grid search, plotting and comparing clusterization parameters)

# 1) best_params_search.py
     - contains a grid which is used in an exhaustive search of parameters. 
     - for each combination of eps and min_samples this code clusterizes data and 
     scores the quality with Silhouette Score and Calinski-Harabazs index. 
     - output file contain a consolidated table with scores and number of clusters
     for corresponding parameters


# 2) Plot_scores.ipynb
     - a note book that reads best_params_search.py outfile (.xlsx) 
     - then heatmaps are plotted. Each cell is colorcoded, and the values inside
     are number of clusters obtained with corresponding parameters pair.


# 3) best_params_clustering_main.py
     - contains both 1 and 2 steps and can be started with in a console


### Plotting a region of a cromosome

To plot a region:

1) You'll need a clusterized chromosome (SNP_clustering_v0.9.1.py)

   python3 SNP_clustering_v0.9.1.py /mnt/wd/nsap/in_data1/chr1_matrix.ld /mnt/wd/nsap/in_data1/chr1.snplist all
   
   written hlist files (hap_file_hdbscan.py -> hlists/)
   python3 hap_file_hdbscan.py
   
   
3) run restore-cluster-membreship.py path_to_hlist path_to_bim path_to_output
   - takes bim file 
   - for every row appends a cluster number if any from hlist file
   - writes outfile (cluster_files/imp1/dbscan or hdbscan)
   
   python3 restore-cluster-membership.py /mnt/wd/nsap/hlists/comp_param_chr1_A_CSM075MS2_hdbscan.hlist /mnt/wd/nsap/imp1/chr1.bim /mnt/wd/nsap/Clustering/plot_region/cluster_files/imp1/hdbscan/chr1.csv


4) run batch-plot-heatmap-01.sh (for imp1) and -02.sh (for imp2)
   bash batch-plot-heatmap-01.sh
