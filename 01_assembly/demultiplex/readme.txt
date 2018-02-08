Additional demultiplexing was done for the TT and GSK fragment libraries since the read count was low.

Specifically: 
for the pooled lane in 141017_D00365_0363_AC57WLANXX -> the script from RC informatics available here (https://github.com/harvardinformatics/bioblog/tree/master/Illumina.Demultiplex) was used
	run time options specified a hamming distance of 3
	these reads that could additionally be recovered are in the demultiplex1 files

for the unpooled lanes on 140925_D00365_0350_AC57WPANXX -> PhiX reads were removed but otherwise everything was retained
	PhiX reads were removed using BBDuk - sbatch scripts in this directory
	the non-phiX reads in the demultiplex2 files
