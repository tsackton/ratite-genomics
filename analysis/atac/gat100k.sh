#!/bin/bash
#SBATCH -J gat
#SBATCH -n 1                     # Use 1 cores for the job
#SBATCH -t 5-00:00                 # Runtime in D-HH:MM
#SBATCH -p general         # Partition to submit to
#SBATCH --mem=20000               # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o gat.%A.out  # File to which STDOUT will be written
#SBATCH -e gat.%A.err  # File to which STDERR will be written

#$1 is name for file
#$2 is space delimited list of inputs

source activate python2_pgrayson

gat-run.py --ignore-segment-tracks --segments=${1} --annotations=sorted_gG4_allTests_annotations_updateBayes.bed --workspace=workspace_galGal4.bed --num-samples=100000 --output-counts-pattern=count_100k_${1}_%s.txt --log=gat_100k_${1}.log > ${1}_100k_gat.out
