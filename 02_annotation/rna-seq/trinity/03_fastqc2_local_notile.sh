#!/bin/bash

#get work dir
WORKDIR=${1:?} #directory to process, exit if not set

#get files
for FQ in $(ls $WORKDIR/trimmed/dronov_sra.R1.pair.fastq.gz);
do
	mv $FQ inseq.fq.gz
	seqtk seq -C inseq.fq.gz > out.fq
	gzip out.fq 
	mv out.fq.gz $FQ
	fastqc -o $WORKDIR/fastqc_reports $FQ
done

