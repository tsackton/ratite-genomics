#!/bin/bash

SPEC=$1

java -jar /n/home12/tsackton/sw/progs/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads	24 orig/${SPEC}_R1_ho_220_1P.fastq.gz orig/${SPEC}_R1_ho_220_2P.fastq.gz \
        reapr/${SPEC}_frag_1.fastq.gz /dev/null	reapr/${SPEC}_frag_2.fastq.gz /dev/null	MAXINFO:100:0.9	CROP:100 MINLEN:100
