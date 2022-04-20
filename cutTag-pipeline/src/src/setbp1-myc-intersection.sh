#!/bin/bash

mkdir -p src/custom

# download files ----------------------------------------------------------------------------------

# wget SRX2085893 from ChIP-Atlas. This is SETBP1 ChIP-Seq from HEK 293 cells aligned to hg38
wget -P src/custom http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/eachData/bed05/SRX2085893.05.bed

# wget ENCFF152JNC from ENCODE. This is MYC ChIP-Seq aligned to mm10.
wget -P src/custom https://www.encodeproject.org/files/ENCFF152JNC/@@download/ENCFF152JNC.bed.gz
zcat src/custom/ENCFF152JNC.bed.gz | sort k1,1 -k2,2n > src/custom/ENCFF152JNC.sorted.bed

# liftOver hg38 --> mm10 chain file
wget -P src/custom https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz

# liftOver setbp1 ---------------------------------------------------------------------------------

liftOver -bedPlus=6 src/custom/SRX2085893.05.bed src/custom/hg38ToMm10.over.chain.gz \
    src/custom/setbp1.hg38.mapped.bed src/custom/setbp1.hg38.unmapped.bed

liftOver -bedPlus=6 src/custom/SETBP1_WT_hg19.bed src/custom/hg19ToMm10.over.chain.gz \
    src/custom/setbp1.hg38.mapped.bed src/custom/setbp1.hg38.unmapped.bed
