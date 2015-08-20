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
generated chicken-trained SNAP models ourselves. See maker/snap/train_snap_chicken.sh for details on
training.

The control files for the MAKER runs are in maker/run1/.

Retraining
----------

We retrained both SNAP and AUGUSTUS based on high-confidence gene models from the initial maker run for
notPer and droNov (the two species we have RNA-seq data for). See the maker/snap/ and maker/augustus/ subdirectories
for code to do this retraining.

Final MAKER run.
----------------

For our final MAKER run, we used the trained gene predictors, as well as expanded evidence sources (TopHat junctions
in additional to Trinity assemblies). The control files for the final MAKER runs are in maker/run2/.


