#!/bin/bash

if [ ! -d log ]; then
  mkdir log
fi

cr="/home/groups/MaxsonLab/kongg/carratt4/cite_seq_primary/cellranger-5.0.1/cellranger"

srun -c 20 --mem=100G --time=24:00:00 $cr count \
  --id=dmso \
  --transcriptome=/home/groups/MaxsonLab/indices/GRch38/refdata-gex-GRCh38-2020-A \
  --libraries=config/dmso.csv \
  --feature-ref=config/feature_ref.csv \
  --localcores=20 \
  --localmem=100 \
  --chemistry=SC3Pv3 > log/dmso.out 2> log/dmso.err &

srun -c 20 --mem=100G --time=24:00:00 $cr count \
  --id=ory \
  --transcriptome=/home/groups/MaxsonLab/indices/GRch38/refdata-gex-GRCh38-2020-A \
  --libraries=config/ory.csv \
  --feature-ref=config/feature_ref.csv \
  --localcores=20 \
  --localmem=100 \
  --chemistry=SC3Pv3 > log/ory.out 2> log/ory.err &
