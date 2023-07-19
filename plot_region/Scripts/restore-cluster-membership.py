#!/usr/bin/env python3
## Restore cluster membership for each SNP from hlist files
# python3 restore-cluster-membership.py /mnt/wd/nsap/Clustering/main_pipeline/hlists/A_chr10_CSM035MS2_EPS075MS2_dbscan.hlist /mnt/wd/nsap/imp1/chr10.bim /mnt/wd/nsap/Clustering/plot_region/cluster_files/imp1/dbscan/A_chr10_CSM035MS2_EPS075MS2_dbscan.csv


import sys
import os.path
import pandas as pd

# Set functions
def proc_bim(x, y):
	df = pd.DataFrame(columns = ('chr', 'pos', 'rs', 'cl'))
	f = open(x, 'r')
	for line in f:
		l = line.strip()
		a = l.split()
		d = {'chr':a[0], 'pos':a[3], 'rs':a[1], 'cl':None}
		if d['rs'] in y:
			d['cl'] = y[d['rs']]
		df = df.append(d, ignore_index = True)
		
	f.close()
	return(df)

# Initiate
if len(sys.argv) < 4:
	print('Input wasn\'t defined')
	sys.exit(1)

hl = sys.argv[1]
bim = sys.argv[2]
output = sys.argv[3]

# Check input
if not os.path.exists(hl):
	print(f'{hl} wasn\'t found')
	sys.exit(1)
if not os.path.exists(bim):
	print(f'{bim} wasn\'t found')
	sys.exit(1)

# Show input
print("Restore cluster membership") 
print("hlist: ", hl)
print("bim: ", bim)
print("output: ", output)

# Initiate empty dictionary for cluster memberships
cl = dict()

# Read hlist file line by line and extract list of SNPs
# belonging to the same cluster
f1 = open(hl, 'r')
for line in f1:
	l1 = line.strip()
	l2 = l1.split()
	l3 = l2[1].split("-")

	# Add new element to dictionary
	k = l3[0]
	if k not in cl:
		cl[k] = l2[2:] 

f1.close()

# Make dictionary of SNPs with cluster memberships
rs = dict()
n = 1
for k in cl:
	for v in cl[k]:
		if v not in rs:
			rs[v] = n
	n += 1

# Process bim file
out = proc_bim(bim, rs)

# Save to file
out.to_csv(output, na_rep='None', index = False)

