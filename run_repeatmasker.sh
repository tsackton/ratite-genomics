#!/bin/bash

GENOME=$1

RepeatMasker -pa 16 -s -norna -species aves -dir ${GENOME}_rm -xsmall -gff $GENOME.fa.gz