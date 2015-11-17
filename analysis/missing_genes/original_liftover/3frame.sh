#!/bin/sh
for FILE in $(ls *.fa)
do
        echo "Submitting job for $FILE"
        sbatch Stopsplitter.sbatch $FILE
        sleep 1
done
