#!/bin/bash

#uses trinity 2.0.6, local install

#get species to process
SP=$1
LAUNCHPATH=$(pwd -P)

# get left and right reads
LEFT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R1.*.fastq.gz | perl -p -i -e 's/\s+//')
RIGHT=$(ls -m $LAUNCHPATH/$SP/trimmed/*.R2.*.fastq.gz | perl -p -i -e 's/\s+//')

# run Trinity
Trinity \
 --seqType fq \
 --max_memory 500G \
 --left $LEFT \
 --right $RIGHT \
 --normalize_reads \
 --output $SP/trinity \
 --CPU 8 \
 --verbose \
 --min_kmer_cov 1 \
 --grid_conf trinity_SLURM_conf.txt \
 --group_pairs_distance 800

