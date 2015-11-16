#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 0-8:00:00 #2-0:00:00 Runtime in days-hours:min:sec
#SBATCH --mem 1200
#SBATCH -p serial_requeue #general #Partition to submit to
#SBATCH -o hlo%A.out
#SBTACH -e hlo%A.err

halLiftover --outPSLWithName --tab /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal galGal Chicken_CDS_HLO.bed $1 $1.psl
