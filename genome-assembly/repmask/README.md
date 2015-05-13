Repeat Masking
===============

We used RepeatMasker version 4.0.5 to mask repetitive elements in our newly sequenced genomes, using NCBI/RMBLAST 2.2.27+ as the search engine and version 20140131 of the RepBase RepeatMasker library.

We called RepeatMasker as follows:
  RepeatMasker -par 8 -species <GROUP> -dir . -gff -excln </PATH/TO/GENOME>

Where <GROUP> is the species library to use, in this case aves, and </PATH/TO/GENOME> points to the genome sequence to mask.