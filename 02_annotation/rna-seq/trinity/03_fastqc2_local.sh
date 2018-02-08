#!/bin/bash

#get work dir
WORKDIR=${1:?} #directory to process, exit if not set

#get files
for FQ in $(ls $WORKDIR/trimmed/*.fastq.gz);
do
	fastqc -o $WORKDIR/fastqc_reports $FQ
done

