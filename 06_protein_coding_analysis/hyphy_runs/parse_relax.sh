

#!/bin/bash
DATE=$(date +%F)
TESTCLASS=$1
WORKDIR="RELAX_"$TESTCLASS
find $WORKDIR -type f -name output -exec grep "Likelihood ratio test" {} \+ >> relax_parsed_Pval_${TESTCLASS}_${DATE} & 
find $WORKDIR -type f -name output -exec grep "Relaxation parameter K" {} \+ >> relax_parsed_K_${TESTCLASS}_${DATE} &
