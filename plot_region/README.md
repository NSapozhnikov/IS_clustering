This folder contains scripts for plotting a region of a cromosome.

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
