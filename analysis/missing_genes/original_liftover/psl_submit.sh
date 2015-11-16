#!/bin/sh
for FILE in $(ls *.fa)
do
	echo "Submitting job for $FILE"
	sbatch psl_dict.sbatch $FILE
	sleep 1
done
