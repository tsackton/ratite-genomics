#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-23:00:00 #2-0:00:00 Runtime in days-hours:min:sec
#SBATCH --mem 1200
#SBATCH -p serial_requeue #general #Partition to submit to
#SBATCH -o hlo_cnee%A.out
#SBTACH -e hlo_cnee%A.err

halLiftover --outPSLWithName --tab /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal galGal most_conserved_final.tree2.bed $1 $1.psl
