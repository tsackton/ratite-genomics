RNA-seq
===============

We used a combination of existing RNA-seq resources from SRA and new RNA-seq data to generate
Trinity assemblies and TopHat alignments for use in MAKER annotation.

*Note: need to add list of SRA accesions*

For each species, we first trimmed and preprocessed sequences (code in trinity subdir: 00_extract_sra.sh,
01_concat.sh, 02_trim.sh, 03_fastqc2.sh), and then ran Trinity to generate de novo assemblies (trinity/04_run_trinity.sh) and
TopHat to map reads (tophat/).
