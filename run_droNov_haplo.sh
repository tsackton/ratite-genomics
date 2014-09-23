#!/bin/bash

#prep path
source setup_allpaths.sh
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=droNov DATA_SUBDIR=RAWSEQ RUN=20140901 HAPLOIDIFY=TRUE THREADS=32 OVERWRITE=TRUE 1> droNov_20140901.log

