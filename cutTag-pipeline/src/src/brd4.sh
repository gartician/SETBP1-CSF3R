#!/bin/bash

# intersect DE peaks with brd4 binding sites from myeloerythroid (MEL) cell line from MOUSE ENCODE.

# just wget ENCFF152JNC and then bedtools sort the file.

mkdir data/brd4

brd4_binding_sites="src/GSM3018457_BRD4_SenDMSO_peaks.narrowPeak"

for file in $(find data/deseq2 -name "*differential*.bed" | sort); do

	contrast=$(basename $file | cut -d- -f1-2)
	mark=$(basename $file | cut -d- -f3)
	dir=$(basename $file | cut -d- -f5)
	sig=$(basename $file | cut -d- -f6 | cut -d. -f1)

	outfile="data/brd4/${contrast}-${mark}-${dir}-${sig}-brd4-intersection.bed"
	bedtools intersect -u -a $file -b $brd4_binding_sites > $outfile

done

conda activate homer

genome="/home/groups/MaxsonLab/indices/mm10/mm10.fa"
sbatch --cpus-per-task=8 -o data/brd4/homer-brd4-H3K27ac-down-05/log --wrap="findMotifsGenome.pl data/brd4/DoxB-DoxWdB-H3K27Ac-down-05-brd4-intersection.bed $genome data/brd4/homer-brd4-H3K27ac-down-05 -size 200 -p 8" &
sbatch --cpus-per-task=8 -o data/brd4/homer-brd4-H3K27ac-up-05/log --wrap="findMotifsGenome.pl data/brd4/DoxB-DoxWdB-H3K27Ac-up-05-brd4-intersection.bed $genome data/brd4/homer-brd4-H3K27ac-up-05 -size 200 -p 8" &