Genome annotation
===============

To annotate protein-coding genes in our newly assembled genomes, we used [MAKER v2.31.8](http://www.yandell-lab.org/software/maker.html).
We used a 2-step approach where we first annotated the emu and Chilean tinamou genomes, 
for which we have RNA-seq data, used the initial annotations to build improved SNAP and AUGUSTUS 
models, and then annotated all 10 species using the improved models. For the species without
same-species RNA-seq data, we use the closest appropriate alternate-species RNA-seq as 
evidence, but we do not retrain gene models for these species. 

RNA-seq data
------------

All code to process RNA-seq data is included in this subdirectory. This includes mapping and Trinity assemblies.


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

IMPORTANT NOTE ADDED 10/22/2015:
--------
We have preliminary evidence that one of the emu tissue samples (emu female liver) is actually a tinamou sample. Implications:
-the Trinity assembly of the droNov reads will contain some tinamou contamination
-the TopHat mapping of the droNov reads will contain some tinamou contamination in the junctions.bed files
-thus, the Trinity assemblies and TopHat mappings of the pooled data should not be used for downstream analysis

However, we believe that the gene models should be relatively robust to a low level of tinamou sequence in the input droNov sequence:
-because MAKER relies on evidence that maps to the genome, divergent transcripts assembled by Trinity that are derived from tinamou sequence generally won't map
-in many cases in the Trinity assembly, the tinamou sequence ends up as an alternative splice form; however, when these are mapped to the genome the "tinamou" isoform either 1) won't map, or 2) will map but overlapping with an "emu" isoform
-finally, we filtered the junctions.bed files from TopHat before using, so junctions with weak support will not be passed to MAKER

Ultimately we see little evidence that the droNov gene models are notably different in any way than the other ratite gene models, and given the computational costs to rerun MAKER we are considering these models frozen, despite being produced with pooled RNA-seq data that included some misclassified sequence (i.e, tinamou reads called as emu)


NCBI
-----

After NCBI submission screening, we generated mask files to screen for gene models that may potentially contain masked sequence, and to update GFF files. Ultimately we did not submit gene models to NCBI but retain the code used to generate updated GFFs in case it is useful. See final_ncbi_submission_code.md for details.


