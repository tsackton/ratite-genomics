#!/bin/bash

makeblastdb -dbtype nucl -in GGALL_genome.fas -parse_seqids
tblastn -db GGALL_genome.fas -query BUSCO_GGALL.fas -evalue 0.1 -out GGALL.blastout.txt -outfmt 6 -seg no -num_threads 24
makeblastdb -in BUSCO_GGALL.fas -parse_seqids 
blastp -db BUSCO_GGALL.fas -query BUSCO_GGALL.fas -evalue 0.1 -out BUSCO_GGALL.SELF.blastout.txt -outfmt 6 -seg no -num_threads 24
perl BUSCO.pl --refine_reference GGALL 95
makeblastdb -in BUSCO_GGALL_refined.fas -parse_seqids 
blastp -db BUSCO_GGALL_refined.fas -query BUSCO_GGALL_refined.fas -evalue 0.1 -out BUSCO_GGALL_refined.SELF.blastout.txt -outfmt 6 -seg no -num_threads 24
tblastn -db GGALL_genome.fas -query BUSCO_GGALL_refined.fas -evalue 0.1 -out GGALL_refined.blastout.txt -outfmt 6 -seg no -num_threads 24
perl BUSCO.pl --run_refined GGALL GGALL

