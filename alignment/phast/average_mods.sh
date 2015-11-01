#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -t 0-02:00
#SBATCH --mem 1000
#SBATCH -J aveRho

ls mods/*.cons.mod > cons.txt
phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
ls mods/*.noncons.mod > noncons.txt
phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod 

