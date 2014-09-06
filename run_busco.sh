#!/bin/bash

ASSEMBLYPATH=$1
SPEC=$2
cp $ASSEMBLYPATH ${SPEC}_genome.fas
makeblastdb -dbtype nucl -in ${SPEC}_genome.fas -parse_seqids 
tblastn -db ${SPEC}_genome.fas -query BUSCO_GGALL_refined.fas -evalue 0.1 -out ${SPEC}_refined.blastout.txt -outfmt 6 -seg no -num_threads 24
perl BUSCO.pl --run_refined $SPEC DMELA

