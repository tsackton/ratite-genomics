#!/bin/bash

#code to clean up paml files for ancestral reconstruction
rm -r anc
rm -r branchsite
rm -r branchsitenull
rm -r branchnull
rm -r branch
mkdir -p anc
mkdir -p branchsite
mkdir -p branchsitenull
mkdir -p branchnull
mkdir -p branch
for FILE in $(find . -name rst | grep "ancrec"); do TMP=${FILE%/HOG*}; HOG=${TMP##*/}; echo $HOG; cp -a $FILE anc/$HOG.paml.events; done
for FILE in $(find . -type f -name "branchsite.out" | grep "branchsite.ctl"); do TMP=${FILE%/HOG*}; HOG=${TMP##*/}; echo $HOG; cp -a $FILE branchsite/$HOG.brstalt; done
for FILE in $(find . -type f -name "branchsitenull.out" | grep "branchsitenull.ctl"); do TMP=${FILE%/HOG*}; HOG=${TMP##*/}; echo $HOG; cp -a $FILE branchsitenull/$HOG.brstnull; done
for FILE in $(find . -type f -name "ancrec.out" | grep "ancrec"); do TMP=${FILE%/HOG*}; HOG=${TMP##*/}; echo $HOG; cp -a $FILE branchnull/$HOG.brnull; done
for FILE in $(find . -type f -name "branch.out" | grep "branch.ctl"); do TMP=${FILE%/HOG*}; HOG=${TMP##*/}; echo $HOG; cp -a $FILE branch/$HOG.bralt; done
./infer_mut_type_v2.pl anc > paml.muts.run1
./get_branch_info.pl branch \.bralt > branch_alt.run1
./get_branch_info.pl branchnull \.brnull > branch_null.run1
./get_branchsite_info.pl branchsite \.brstalt > branchsite_alt.run1
./get_branchsite_info.pl branchsitenull \.brstnull > branchsite_null.run1
