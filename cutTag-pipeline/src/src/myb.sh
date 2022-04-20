#!/bin/bash

# intersect DE peaks with myb binding sites from myeloerythroid (MEL) cell line from MOUSE ENCODE.

# just wget ENCFF911NHJ and then bedtools sort the file.

mkdir data/myb

myb_binding_sites="src/custom/ENCFF911NHJ.sorted.bed"

for file in $(find data/deseq2 -name "*differential*.bed" | sort); do

        contrast=$(basename $file | cut -d- -f1-2)
        mark=$(basename $file | cut -d- -f3)
        dir=$(basename $file | cut -d- -f5)
        sig=$(basename $file | cut -d- -f6 | cut -d. -f1)

        outfile="data/myb/${contrast}-${mark}-${dir}-${sig}-myb-intersection.bed"
        bedtools intersect -u -a $file -b $myb_binding_sites > $outfile

done
