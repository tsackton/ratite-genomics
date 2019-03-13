#clean up existing modules
module purge
#setup py3 conda environment
module load gcc/7.1.0-fasrc01 
module load Anaconda3/5.0.1-fasrc02
source activate py3
#load modules, guarantees versions will persist as these module versions
#if you want to use conda versions need to comment out
module load jdk/1.8.0_45-fasrc01
module load samtools/1.5-fasrc02
module load seqtk/1.2-fasrc01
module load bedtools2/2.26.0-fasrc01
module load ucsc/20150820-fasrc01
module load R/3.5.0-fasrc02
#fix R LIBS
export R_LIBS_USER=~/sw/R/3.5/:$R_LIBS_USER
#fix LD LIBRARY PATH for clusterProfiler
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/n/home12/tsackton/.conda/envs/py3/lib
#get rid of 'py3' prepended
export PS1='[\u@\h \W]\$ '


