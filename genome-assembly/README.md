Genome assembly
===============

- Tim Sackton (tsackton at oeb.harvard.edu)
- Michele Clamp (mclamp at g.harvard.edu)
- Scott Edwards (sedwards at oeb.harvard.edu)

1. Prior to assembly, we preprocessed the sequencing reads using trim_fastq.sh
2. After preprocessing, we verified quality by generating FastQC reports with run_fastqc.sh
3. We assembled genomes with ALLPATHS-LG using final_assembly.sh

Subdirectories:

- logs: contains ALLPATHS-LG output logs for each species we assembled
- libraries: contains ALLPATHS-LG libs.csv and groups.csv files for each species
- qc: contains scripts and other materials for running QC on the assemblies