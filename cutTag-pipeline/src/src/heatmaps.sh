#!/bin/bash

sig="05" # 05 or 01
banlist=$1
# banlist="/home/groups/MaxsonLab/indices/mm10/mm10.blacklist.v2.bed"
mkdir data/computeMatrix

all_marks=$(cut -f4 samplesheet.tsv | tail -n+2 | grep -v "IgG" | sort | uniq)

# assume merge_bw rule is finished
for mark in $all_marks; do

	# identify your sample
	bigwigs=$(find data/mergebw -name "*${mark}*" | tr '\n' ' ')
	up="data/deseq2/${mark}/*${mark}*-up-${sig}.bed"
	down="data/deseq2/${mark}/*${mark}*-down-${sig}.bed"
	consensus="data/counts/${mark}_consensus.bed"

	# file I/O
	consensus_output="data/computeMatrix/${mark}.consensus.gz"
	differential_output="data/computeMatrix/${mark}.differential.gz"

	# consensus signal
	echo "Making matrix for $bigwigs"
	echo "computeMatrix reference-point --referencePoint center -p 16 -S $bigwigs -R $consensus -a 3000 -b 3000 -o $matrix -bl $banlist --smartLabels"
	sbatch -c 16 --mem 8G --wait --wrap "computeMatrix reference-point --referencePoint center -p 16 -S $bigwigs -R $consensus -a 3000 -b 3000 -o $consensus_output -bl $banlist --smartLabels" &

	# differential signal
	echo "computeMatrix reference-point --referencePoint center -p 16 -S $bigwigs -R $consensus -a 3000 -b 3000 -o $matrix -bl $banlist --smartLabels"
	sbatch -c 16 --mem 8G --wait --wrap "computeMatrix reference-point --referencePoint center -p 16 -S $bigwigs -R $up $down -a 3000 -b 3000 -o $differential_output -bl $banlist --smartLabels" &

done
wait

mkdir -p data/figures/heatmaps
for matrix in $(find data/computeMatrix -name "*.gz"); do

	# identify your sample
	mark_type=$(basename $matrix | cut -d. -f1-2) # "H3K27Ac.consensus" or "H3K27Ac.differential"

	# file I/O
	outfile="data/figures/heatmaps/${mark_type}.pdf"

	# make the heatmaps
	plotHeatmap -m $matrix -o $outfile &

done
wait