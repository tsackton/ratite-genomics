#!/bin/bash

#prep path
source setup_allpaths.sh
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=aptHaa DATA_SUBDIR=data RUN=20141028 HAPLOIDIFY=TRUE THREADS=32 OVERWRITE=TRUE > aptHaa_20141028.log

