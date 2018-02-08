#!/bin/bash
#SBATCH -p general
#SBATCH --mem-per-cpu 6000
#SBATCH -n 128
#SBATCH -t 3-00:00
#SBATCH -o notper_mpi5.out
#SBATCH -e notper_mpi5.err
#SBATCH -J notperMPI

source new-modules.sh
module load openmpi/1.8.3-fasrc02
export LD_PRELOAD=/n/sw/fasrcsw/apps/Comp/gcc/4.8.2-fasrc01/openmpi/1.8.3-fasrc02/lib64/libmpi.so
mpiexec  -mca btl ^openib -np 128 maker -fix_nucleotides
