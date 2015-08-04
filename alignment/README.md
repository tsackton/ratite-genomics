Whole genome alignment
===============

We used [progressiveCactus](https://github.com/glennhickey/progressiveCactus) to generate
a whole genome alignment from 35 bird genomes and 7 non-avian reptile outgroups. Our
version of progressiveCactus uses a modified version of [jobTree](https://github.com/harvardinformatics/jobTree) in order to work with SLURM.

Alignment steps
---------------

1. Selection of genomes: in addition to our 10 new palaeognath genomes, we focused on
other bird genomes of relatively high quality measured by both contiguity statistics 
(contig and scaffold N50) and completeness (BUSCO). We make a few exceptions for important 
outgroups. The full list of species used, including assembly versions, is in 
species_list.tsv, and the script to fetch them from NCBI is get_genomes.sh.

2. Guide tree generation: see branchlengths subdirectory

3. Whole-genome alignment: see progressiveCactus subdirectory
