#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-6:00:00 #0-0:05:00 Runtime in days-hours:min:sec
#SBATCH --mem 2048
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o hal2maf%A_%a.out
#SBTACH -e hal2maf%A_%a.err

mkdir -p output/${SLURM_ARRAY_TASK_ID}
mkdir -p log/${SLURM_ARRAY_TASK_ID}
cd Beds/${SLURM_ARRAY_TASK_ID}
for FILE in $(ls *.bed) 
do
        echo ${FILE}_look
        if [ ! -f ../../log/${SLURM_ARRAY_TASK_ID}/${FILE}.done ]
        then
                echo ${FILE}_notFound
                echo ${FILE}_run
                REFCHR=$(head -n 1 ${FILE} | cut -f1,1) #necessary because hal2maf requires refSequence
                hal2maf --refGenome galGal --refTargets ${FILE} --refSequence ${REFCHR} --maxRefGap 50000 --noAncestors /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal ../../output/${SLURM_ARRAY_TASK_ID}/${FILE}.maf
                touch ../../log/${SLURM_ARRAY_TASK_ID}/${FILE}.done
                echo ${FILE}_touched
        fi
done
