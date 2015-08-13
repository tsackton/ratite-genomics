MAKER run2
---------

For the second (final) run of MAKER, we use the following inputs:
1. RNA-seq data: TopHat junctions converted to GFF (see TopHat directory)
2. RNA-seq data: Trinity assemblies (see Trinity directory)
3. Protein data: full proteomes from 10 species (same as MAKER run 1)
4. Trained gene finders: SNAP & Augustus (notPer training for tinamous, droNov training for ratites)


Setup:
base directory: /n/regal/edwards_lab/ratites/maker2
rnaseq: $BASE/evidence/rnaseq (includes GFF files and Trinity assemblies)
snap hmms: $BASE/evidence/snap
protein evidence: $BASE/evidence/prot
genomes: $BASE/inputs
run directories: $BASE/annotation/$SPECIES

setup_maker.sh sets up directories (but still has paths hardcoded for copying files)
