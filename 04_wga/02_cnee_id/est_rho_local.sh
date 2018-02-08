#!/bin/bash

for FILE in $(cat files)
do
        SAMP=${FILE%.ss}
        STATUS=$(tail -n 1 ./logs/$SAMP.out)
        if [[ "$STATUS" != "Done." ]]
        then
		phastCons --expected-length=45 --target-coverage=0.3 --estimate-rho ./mods/$SAMP.rho --no-post-probs --msa-format SS --log ./logs/$SAMP.log ../chunks/$SAMP.ss ../../../neutMods/neut_ver1_final.named.mod &> ./logs/$SAMP.out &
	fi
done
