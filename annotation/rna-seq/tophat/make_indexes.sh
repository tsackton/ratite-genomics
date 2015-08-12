#!/bin/bash

cp ~/ratite_store/genomes/ASSEMBLY_FREEZE/*.fa.gz .
gunzip *.gz
module load tophat
for GENOME in $(ls *.fa)
do
	 bowtie2-build $GENOME ${GENOME%%.*} &
done

