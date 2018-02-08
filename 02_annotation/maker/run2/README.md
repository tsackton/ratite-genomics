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

RNA-seq evidence logic
----------------------

For each species, we attempt to use the RNA-seq evidence that is phylogenetically closest, as follows:

1. Kiwis: we treat the aptMan data (Trinity assembly + TopHat junctions) as same-species data, and use 
droNov data for the alt-ests
2. Emu: the droNov data is same species, we use kiwi data as alt-ests
3. Cassowary: we do not have same species data, and since both kiwi and emu are relatively closely related to cassowary 
(in the same clade) we use both kiwi and emu data for alt-ests.
4. Rheas: rheas are tricky, as a) their phylogenetic placement is uncertain and b) we have no closely related RNA-seq evidence. 
So in this case, we used the droNov/kiwi merged data for the TopHat mapping, and the droNov data only for the alt-est. But we acknowledge
this choice is a bit arbitrary; however, the computational cost of MAKER is too high to repeat with alternative choices to see what
impact it has.
5. Tinamous: we used the tinamou data only (notPer), either as same-species or alt-species as appropriate.

To submit the jobs, we ran:
```for SP in aptHaa aptRow aptOwe droNov casCas cryCin notPer eudEle rheAme rhePen; do sbatch maker_mpi.sh $SP; done```
from the MAKER base directory

Cleanup
-------

After running MAKER, we clean up and parse the results with this code:
 ```while read LINE; do ./cleanup_maker.sh $LINE & done < maker_key.txt```
 