#!/usr/bin/env python3
# Plot LD, results of clusterization and genome
# Pillow documentation
# https://pillow.readthedocs.io/en/stable/reference/ImageDraw.html
# Author: Gennady Khvorykh, info@inzilico.com
# Date: August 25, 2022

import math
import argparse
import sys
import os.path
import re
import time
import functools as ft
import pandas as pd
import numpy as np
from PIL import Image, ImageDraw, ImageFont
import matplotlib.pyplot as plt

## Start functions

def check_region(x):
	try:
		chrom, start, end = re.split(':|-', x)
		start = int(start)
		end = int(end)
	except:
		print(f'Wrong region: {x}')
		sys.exit(1)
	if 'chr' not in chrom:
		print(f'Wrong region: {chrom} doesn\'t have \'chr\' prefix')
		sys.exit(1)
	if start >= end:
		print(f'Wrong region: {start} isn\'t less {end}')
		sys.exit(1)
	return([chrom, int(start), int(end)])

### Original
def plot_heatmap(m, L, W, a, margins, draw, ld_show = False):
	# Initiate color pallete
	pal = [(0, 85, 63),(0, 132, 98),(0, 180, 132),(0, 227, 167),(19, 255, 193),
				(66, 255, 205),(113, 255, 218),(161, 255, 230),(208, 255, 243),(255, 255, 255)]
	palr = list(reversed(pal))

	# Get the number of markers
	n = m.shape[0]

	# Make ticks
	#ticks = []
	#h = 10
	#for i in range(n):
	#		x = 50+i*2*a
	#		y = W-50
	#		ticks.append([x, y, x, y+h])

	# Plot line
	#draw.line((50, W-50, L-55, W-50), fill=(0, 0, 0), width=2)

	#for i in range(n):
	#		draw.line(ticks[i], fill=(0, 0, 0), width=1)
	#		draw.text((ticks[i][2]-5, ticks[i][3]+2), str(ind[i]), font = font, fill=(0, 0, 0))

	# Add square diamonds
	s = 1
	for i in range(n-1):
			l = 0
			for j in range(i, n-1):
					x = margins['left']+a+2*j*a-i*a
					y = a*n+margins['top']-i*a
					# Add square diamonds
					v1 = m.iloc[l, l+s]
					ic = math.floor(v1*10)
					if ic >= 10:
						ic = 9
					draw.regular_polygon((x, y, a), n_sides=4, rotation=45, fill=palr[ic], outline=(0, 85, 63))
					if ld_show:
						v2 = round(v1, 1)
						draw.text((x-5, y-4), str(v2), fill=(0, 0, 0))
					l = l + 1
			s = s + 1

	return draw

### With the colorbar does not work properly
# def plot_heatmap(m, L, W, a, margins, draw, ld_show = False):
#     # Initiate color pallete
#     pal = [(0, 85, 63),(0, 132, 98),(0, 180, 132),(0, 227, 167),(19, 255, 193),
#            (66, 255, 205),(113, 255, 218),(161, 255, 230),(208, 255, 243),(255, 255, 255)]
#     palr = list(reversed(pal))
# 
#     # Get the number of markers
#     n = m.shape[0]
#     # Add square diamonds
#     s = 1
#     for i in range(n-1):
#         l = 0
#         for j in range(i, n-1):
#             x = margins['left']+a+2*j*a-i*a
#             y = a*n+margins['top']-i*a
#             # Add square diamonds
#             v1 = m.iloc[l, l+s]
#             ic = math.floor(v1*10)
#             if ic >= 10: 
#                 ic = 9 
#             draw.regular_polygon((x, y, a), n_sides=4, rotation=45, fill=palr[ic], outline=(0, 85, 63))
#             if ld_show:
#                 v2 = round(v1, 1)
#                 draw.text((x-5, y-4), str(v2), fill=(0, 0, 0))
#             l = l + 1
#         s = s + 1
#     # Add colorbar
#     x = margins['left'] + W + a
#     y0 = margins['top'] + L/2
#     y1 = margins['top'] + L/2 + a*n/2
#     for i, c in enumerate(palr):
#         y = y0 + i*(y1-y0)/10
#         draw.rectangle((x, y, x+a, y+(y1-y0)/10), fill=c)
#     
#     return draw


def add_clusters(m, cl, W, a, margins, s, l, draw):
	# Initiate
	dc = [(230, 25, 75), (60, 180, 75), (255, 225, 25), (0, 130, 200), (245, 130, 48), (145, 30, 180), 
				(70, 240, 240), (240, 50, 230), (210, 245, 60), (250, 190, 212), (0, 128, 128), (220, 190, 255), 
				(170, 110, 40), (255, 250, 200), (128, 0, 0), (170, 255, 195), (128, 128, 0), (255, 215, 180), 
				(0, 0, 128), (128, 128, 128), (255, 255, 255), (0, 0, 0)]
	
	font = ImageFont.truetype("Ubuntu-R.ttf", 15)

	# Set distinct colors for clusters	
	un = cl.unique()
	k = 0
	c = dict()
	for i in un:
		if i == 'None':
			c[i] = 20
		else: 
			if k == 20:
				print("More than 20 clusters")
				sys.exit(1) 
			c[i] = k
			k += 1

	# Add cluster memberships
	r = math.sqrt(2)*a # Radius to plot squares
	n = m.shape[0]
	for i in range(n):
		x = margins['left']+2*a*i
		y = a*n+margins['top']+2*a+s*2*a  
		key = cl.iloc[i]
		fill = dc[c[key]] 
		draw.regular_polygon((x, y, r), n_sides=4, fill=fill, outline=(0, 85, 63))
		if i == (n-1):
			draw.text((x+2*a, y-a/2), l, font = font, fill=(0, 0, 0))
	
	return draw 

def add_labels(m, W, a, margins, s, draw):
	font = ImageFont.truetype("Ubuntu-R.ttf", 15)
	# Get index for lables
	ind = m.index
	n = m.shape[0] 
	#	Add marker labels
	for i in range(n):
		x = margins['left']+2*a*i-5 
		y = a*n+margins['top']+2*a+s*2*a  
		draw.text((x, y), str(ind[i]), font = font, fill=(0, 0, 0))
	return draw


def merge_clusters(clusters):
	dfs = []
	nl = []
	print('Label\tPath\tLines\tClusters\tNone, %')
	for x in clusters:
		l, p = x.split(':')
		if not os.path.isfile(p):
			print(f'{p} doesn\'t exist')
			sys.exit(1)
		df = pd.read_csv(p)
		n = df.shape[0]
		nl.append(n)
		na = round(sum(df['cl'] == 'None')/n*100, 1)
		un = len(df['cl'].unique()) - 1 
		print(f'{l}\t{p}\t{n}\t{un}\t{na}')
		df.rename(columns = {'cl': l}, inplace = True)
		dfs.append(df)
	# Merge all dataframes
	out = ft.reduce(lambda left, right: pd.merge(left, right, on = ('chr', 'pos', 'rs'), how = 'inner'), dfs)
	n = out.shape[0]
	nl.append(n)
	print(f'Merged dataframe: {n} rows')
	res = all(x == nl[0] for x in nl)
	if not res:
		print("Different number of rows!")
		sys.exit(1)
	return out

def load_ld(ld, folder):
	filename = os.path.basename(ld)
	h5 = f'{folder}/{filename}.h5'
	if os.path.isfile(h5): 
		print(f'Loading {h5}')
		r2 = pd.read_hdf(h5)		
	else:
		print(f'Loading {ld}')
		r2 = pd.read_csv(ld, sep = "\t", header = None)
		print(f'Writing {h5}')
		r2.to_hdf(h5, key = 'r2', mode = 'w')
	return r2

def main():
	# Get time
	t1 = time.time()

	# Parse command line arguments
	parser = argparse.ArgumentParser(description='Visulize the results of clusterization')
	parser.add_argument('-r', '--region', help = 'region chrN:start-end', required = True)
	parser.add_argument('-l', '--ld', help = 'path/to/filename.txt with LD matrix', required = True)
	parser.add_argument('-c', '--cluster', help = 'path/to/file.csv with cluster memberships', 
											action = 'extend', nargs = '+', required = True)
	parser.add_argument('-o', '--out', help = 'path/to/output folder to save graphics created as region.png', required = True)
	args = parser.parse_args()

	# Initiate variables
	h5_folder = '/mnt/wd/nsap/imp1/hdf5'
	region = args.region
	ld = args.ld
	clusters = args.cluster
	out = args.out

	# Check input
	chrom, start, end = check_region(region)
	print(f'Region: {chrom}:{start}-{end}')
	cl1 =  merge_clusters(clusters)

	if not os.path.isfile(ld):
		print(f'{ld} doesn\'t exist')
		sys.exit(1)

	out = os.path.normpath(out)
	if not os.path.isdir(out):
		print(f'{out} to save output doesn\'t exist')
		sys.exit(1)

	# Load LD matrix
	r2 = load_ld(ld, h5_folder)
	
	# Subset region by positioins
	f1 = cl1['pos'] >= start 
	f2 = cl1['pos'] <= end
	cl2 = cl1[f1 & f2]
	ind1, ind2 = cl2.index[0], cl2.index[-1] + 1

	# Subset a region from LD matrix by indexes
	m1 = r2.iloc[ind1:ind2, ind1:ind2]

	# Set image size
	n = m1.shape[0]
	a = 15 # Radius to plot diamonds
	margins = {'left':50, 'right':100, 'top':50, 'bottom':50}
	nc = len(clusters)+1 # The number of bottom tracks including labels 
	L = 2*a*n+margins['left']+margins['right']
	W = a*n+margins['top']+margins['bottom']+2*a*nc

	print(f'n: {n}\nL: {L}\nW: {W}')
	
	# Initiate image
	im = Image.new('RGB', (L, W), (255, 255, 255))
	draw = ImageDraw.Draw(im)
	
	# Plot heatmap
	draw = plot_heatmap(m1, L, W, a, margins, draw)

	# Add cluster memberships
	for i in range(len(clusters)):
		l, p = clusters[i].split(':')
		cl = cl2[l]
		draw = add_clusters(m1, cl, W, a, margins, i, l, draw)

	# Add labels
	draw = add_labels(m1, W, a, margins, nc-1, draw)

	# Save the graphics to file
	fn =	f'{out}/{chrom}-{start}-{end}.png'
	im.save(fn)
	print(f'{fn} was saved.')
	
	# Show time elapsed
	t2 = time.time()
	dur = time.strftime("%H:%M:%S", time.gmtime(t2-t1))
	print("Time elapsed: ", dur)

## End functions

if __name__ == "__main__":
	main()
