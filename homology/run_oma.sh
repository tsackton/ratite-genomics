#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 16000
#SBATCH -J oma_ratite
#SBATCH -e oma_%j.err
#SBATCH -o oma_%j.out
#SBATCH --array=1-1000
#SBATCH -t 3-00:00
#SBATCH -p general
#SBATCH --account rc_admin

/n/home12/tsackton/sw/source/OMA.1.0.0/bin/oma -s
