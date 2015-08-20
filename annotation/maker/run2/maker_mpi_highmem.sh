#!/bin/bash
#SBATCH -p general
#SBATCH --mem-per-cpu 8000
#SBATCH -n 256
#SBATCH -t 6-00:00
#SBATCH -o maker_%A.out
#SBATCH -e maker_%A.err
#SBATCH -J mkerMPI

#set up and log species
SP=$1

echo "Running MAKER on $SP."

#load modules
source new-modules.sh
module load openmpi/1.8.3-fasrc02
export LD_PRELOAD=/n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/openmpi/1.8.3-fasrc02/lib64/libmpi.so

#change to correct directory
cd /n/regal/edwards_lab/ratites/maker2/annotation/$SP

#run MAKER
mpiexec -mca btl ^openib -np 256 maker maker_opts_$SP.ctl maker_bopts.ctl maker_exe.ctl -fix_nucleotides
