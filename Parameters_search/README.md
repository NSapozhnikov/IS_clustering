### This folder contains scripts for grid search, plotting and comparing
    clusterization parameters

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
