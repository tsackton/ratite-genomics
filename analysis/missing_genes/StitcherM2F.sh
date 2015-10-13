#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-10:00:00 #0-5:00:00 Runtime in days-hours:min:sec
#SBATCH --mem=400
#SBATCH -p serial_requeue #Partition to submit to
#SBATCH -o minuso/m2fStitcher%A_%a.out
#SBTACH -e minuse/m2fStitcher%A_%a.err

mkdir -p output_mfa/${SLURM_ARRAY_TASK_ID}
mkdir -p log_mfa/${SLURM_ARRAY_TASK_ID}
mkdir -p temp/${SLURM_ARRAY_TASK_ID}
mkdir -p lists

for FILE in $(ls ../../Final/output/${SLURM_ARRAY_TASK_ID}/) #for every file in folder array number
do
	echo ${FILE}_look
        if [ ! -f log_mfa/${SLURM_ARRAY_TASK_ID}/${FILE}.done ]
        then
                echo ${FILE}_notFound
                echo ${FILE}_run
		export FILE #to be used in python
                export SLURM_ARRAY_TASK_ID  #same as above
		python Stitch_wrapper.py #run script
		mafToFastaStitcher --maf ../../Final/output/${SLURM_ARRAY_TASK_ID}/${FILE} --seqs lists/list${SLURM_ARRAY_TASK_ID}.txt --breakpointPenalty 5 --interstitialSequence 20 --outMfa output_mfa/${SLURM_ARRAY_TASK_ID}/${FILE}.fasta
                touch log_mfa/${SLURM_ARRAY_TASK_ID}/${FILE}.done
                echo ${FILE}_touched
       		rm temp/${SLURM_ARRAY_TASK_ID}/* #clear the associated temp folder
	fi
done

