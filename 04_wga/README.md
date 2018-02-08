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
species_list.tsv, and the script to fetch them from NCBI is get_genomes.sh (in the progressiveCactus subdirectory).

2. Guide tree generation: the code in the branchlengths subdirectory includes the methods and data sources for generating the guide tree
we used for our whole genome alignment. The R script compute_br_for_cactus.R contains all the information and scripts needed to
generate our final branch lengths.

3. Whole-genome alignment: the code in the progressiveCactus subdirectory includes our configuration files and
scripts to run progressiveCactus. For more details about progressiveCactus, please see [here](https://github.com/glennhickey/progressiveCactus).

Note about ostrich
--------------

We used the superscaffolds produced by BGI for our whole genome alignment. All other species are either from NCBI or newly assembled for this project.

CNEEs
--------

We identified CNEEs using PHAST (02_cnee_id). These were then loosely assigned to genes based on the nearest transcription start site (03_cnee_annotation).