### This folder contains scripts for clustering -> hlists preparation ->
    haplotype testing -> consolidationg haplotypes after Bonferroni.
## Author: Nikita Sapozhnikov, nikita.sapozhnikov1@gmail.com
## Date: July 17, 2023

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

# To automatically perform steps 1-5 on all chromosomes use all_chr_pipe.py. But keep in mind
  the clusterization parameters.
     python3 all_chr_pipe.py
