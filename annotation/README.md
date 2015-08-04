Genome annotation
===============

To annotate protein-coding genes in our newly assembled genomes, we used [MAKER v2.31.8](http://www.yandell-lab.org/software/maker.html).
We used a 2-step approach where we first annotated the emu and Chilean tinamou genomes, 
for which we have RNA-seq data, used the initial annotations to build improved SNAP and AUGUSTUS 
models, and then annotated all 10 species using the improved models. For the species without
same-species RNA-seq data, we use the closest appropriate alternate-species RNA-seq as 
evidence, but we do not retrain gene models for these species. 

Initial MAKER run.
------------------

For the initial MAKER run, we used the Augustus chicken models distributed with Augustus, and
generated chicken-trained SNAP models ourselves, as follows:


