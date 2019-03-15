Genome assembly
===============

1. Prior to assembly, we preprocessed the sequencing reads using trim_fastq.sh
2. After preprocessing, we verified quality by generating FastQC reports with run_fastqc.sh
3. We assembled genomes with ALLPATHS-LG using final_assembly.sh

Subdirectories:

- kmer: contains the output of the ALLPATHS-LG kmer analysis
- libraries: contains ALLPATHS-LG libs.csv and groups.csv files for each species
- logs: contains ALLPATHS-LG output logs for each species we assembled
- busco: contains scripts and other materials for running BUSCO on the assemblies, as well as summary outputs
- demultiplex: additional scripts for custom demultiplexing of a handful of libraries
- ncbi: workflow / code for fixing errors in genome assemblies identified by NCBI submission pipeline