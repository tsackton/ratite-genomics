#!/bin/sh
for FILE in $(ls *.fa)
do
        echo "Submitting job for $FILE"
        sbatch Into_t_t_launcher.sbatch $FILE
        sleep 1
done
