#!/bin/bash

# intersect DE peaks with MYC binding sites from myeloerythroid (MEL) cell line from MOUSE ENCODE.

# just wget ENCFF152JNC and then bedtools sort the file.

mkdir data/myc

myc_binding_sites="src/custom/ENCFF152JNC.sorted.bed"

for file in $(find data/deseq2 -name "*differential*.bed" | sort); do

	contrast=$(basename $file | cut -d- -f1-2)
	mark=$(basename $file | cut -d- -f3)
	dir=$(basename $file | cut -d- -f5)
	sig=$(basename $file | cut -d- -f6 | cut -d. -f1)

	outfile="data/myc/${contrast}-${mark}-${dir}-${sig}-myc-intersection.bed"
	bedtools intersect -u -a $file -b $myc_binding_sites > $outfile

done