#!/bin/bash

#prep path
source setup_allpaths.sh
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=aptOwe DATA_SUBDIR=RAWSEQ RUN=20140901 HAPLOIDIFY=TRUE THREADS=32 OVERWRITE=TRUE 1> aptOwe_20140901.log

