#!/bin/bash

#prep path
source setup_allpaths.sh
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=rheAme DATA_SUBDIR=RAWSEQ RUN=20140828 HAPLOIDIFY=TRUE THREADS=32 OVERWRITE=TRUE 1> rheAme_20140828.log

