#!/bin/bash

#prep path
source setup_allpaths.sh
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=droNov DATA_SUBDIR=RAWSEQ RUN=20140820 THREADS=56 1> droNov_20140820.log 2> droNov_20140820.err
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=aptOwe DATA_SUBDIR=RAWSEQ RUN=20140820 THREADS=56 1> aptOwe_20140820.log 2> aptOwe_20140820.err
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=rheAme DATA_SUBDIR=RAWSEQ RUN=20140820 THREADS=56 1> rheAme_20140820.log 2> rheAme_20140820.err
RunAllPathsLG PRE=/scratch/tsackton/allpaths_runs REFERENCE_NAME=casCas DATA_SUBDIR=RAWSEQ RUN=20140820 THREADS=56 1> casCas_20140820.log 2> casCas_20140820.err

