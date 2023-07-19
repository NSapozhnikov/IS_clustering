regions=(
	chr1:49210570-50694740
	chr1:23429122-24453078
)


for region	in "${regions[@]}"; do
	IFS=: read -r k reg <<< "${region}"
	
	./plot-heatmap-02.py \
		--region ${k}:${reg} \
		--ld /mnt/wd/nsap/in_data1/${k}_matrix.ld \
		--cluster dbscan:/mnt/wd/nsap/Clustering/plot_region/cluster_files/imp1/dbscan/A_${k}_CSM01MS2_EPS01MS2_dbscan.csv \
		--cluster hdbscan:/mnt/wd/nsap/Clustering/plot_region/cluster_files/imp1/hdbscan/A_${k}_CSM01MS2_EPS01MS2_hdbscan.csv \
		--out /mnt/wd/nsap/Clustering/plot_region/png/A_chr1_CSM01MS2_EPS01MS2/

	echo
done
