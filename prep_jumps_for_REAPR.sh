#!/bin/bash

SPEC=$1
cat orig/${SPEC}3kb.trim_1P.fastq.gz orig/${SPEC}_R1_ho_3kb_1P.fastq.gz | seqtk seq -r - > ./reapr/${SPEC}_jump_1_forREAPR.fastq &
cat orig/${SPEC}3kb.trim_2P.fastq.gz orig/${SPEC}_R1_ho_3kb_2P.fastq.gz | seqtk seq -r - > ./reapr/${SPEC}_jump_2_forREAPR.fastq &
wait
