ratite-genomics
===============

Scripts and other material for the assembly and analysis of ratite genomes.

Preprocessing
-------------

1. rename_files_*.sh: readname raw sequencing read files for ease of use
2. trim_fastq.sh: trim reads with Trimmomatic to remove adapter sequence
3. run_fastqc.sh: run FastQC for quality checking
4. prep_allpaths_runs.sh: prepare sequencing data for AllPaths-LG 
	(note that this script is run on the high memory server and does not use SLURM)
	
Assembly
--------

1. run_allpaths_ratites.sh: script to run initial assemblies for all species
2. run_*_haplo.sh: script to run haploidify=T versions of assemblies for all species
3. source_allpaths.sh: script to set up Allpaths environment on the high memory server (called by above scripts)

Assembly QC
-----------

1. prep_*.sh: preparation scripts for running REAPR and BUSCO
2. run_reapr.sh: run REAPR pipeline
3. run_busco.sh: run BUSCO pipeline
4. run_cegma.sh: run CEGMA pipeline


