regions=(
	chr22:19149580-19268873
	chr22:19279641-19524902
)


for region	in "${regions[@]}"; do
	IFS=: read -r k reg <<< "${region}"
	
	./plot-heatmap-02.py \
		--region ${k}:${reg} \
		--ld /mnt/wd/nsap/in_data2/${k}_matrix.ld \
		--cluster dbscan:/mnt/wd/nsap/ld2/cluster_files/dbscan/${k}.csv \
		--cluster hdbscan:/mnt/wd/nsap/ld2/cluster_files/hdbscan/${k}.csv \
		--out /mnt/wd/nsap/ld2/png/ 

	echo
done
