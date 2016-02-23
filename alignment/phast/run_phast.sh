#CODE TO RUN PHAST AND RELATED ANALYSIS
#MOSTLY RUN INTERACTIVELY ON A LARGE MEMORY MACHINE, SO FEW SLURM SCRIPTS

### GETTING NEUTRAL MODELS ###

#Step 1. Make tree
halStats --tree /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal > tree1.nh
perl -p -i -e 's/Anc\d+//g' tree1.nh 
nw_topology tree1.nh > ratiteTree.nh
nw_labels ratiteTree.nh > species_list

#resolving bifurcations:
#1: reptiles -- turtles as outgroup to archosaurs, gharial + crocs --> croc genome paper
#2: palaeognaths -- several options in different trees
#3: passerines -- ground tit as outgroup to other passerines per Alison's UCE tree
#4: accept Afroaves to resolve landbird polytomy
#5: accept Columbea as sister to passera
#6: balReg + chaVoc = clade (Gruimorphae)
#7: accept Otidae as outgroup to other Passera
#8: Gruimorphae outgroup to waterbirds + landbirds

#now rheas
#ver1 = UCE tree (rheas + tinamous)
#ver2 = Mitchell tree (rheas outgroup to non-ostrichs)
#ver3 = rheas + ECK clade

#Step 2. Get 4-fold degenerate sites based on galGal4 NCBI annotations
#convert to GTF
module load cufflinks
gffread --no-pseudo -C -T -o galGal4.gtf GCF_000002315.3_Gallus_gallus-4.0_genomic.gff 
grep -v "^NC_001323.1" galGal4.gtf | grep "protein_coding" | grep -P "\tCDS\t" > galGal4_filt.gtf

#convert to genepred
gtfToGenePred galGal4_filt.gtf galGal4.gp

#convert to bed
genePredToBed galGal4.gp galGal4.bed

#extract 4D sites
hal4dExtract --conserved --inMemory /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal galGal galGal4.bed galGal4_4d.bed

#not going to use the wrapper scripts as they seem to do odd things. So let's first get a chicken-referenced MAF
mkdir extract_maf
cd extract_maf
cp ../galGal4_4d.bed .
hal2mafMP.py --numProc 48 --splitBySequence --refGenome galGal --noAncestors --noDupes --refTargets galGal4_4d.bed /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal neut4d_input_galGal_ref_ver09222015.maf

#get scaffold list to check that process works
cut -f1,1 galGal4_4d.bed | sort | uniq > galGal_scaffold_list
ls extract_maf/*.maf | wc -l

#fix MAF with sed 
#run the big scaffolds in parallel
for MAF in $(ls extract_maf/*NC*.maf);
do
	sed -i -e 2d $MAF &
done

#run the small scaffolds in serial
for MAF in $(ls extract_maf/*NW*.maf);
do
	sed -i -e 2d $MAF
done

#merge MAFs and make SS file
SPECIES=$(nw_labels ratiteTree.nh | sort | tr '\n' ',')
MSAFILES=$(ls extract_maf/*.maf)

#NOTE: this step removes soft-masked 4d sites##

msa_view --aggregate ${SPECIES%?} --in-format MAF --out-format SS --unordered-ss $MSAFILES > neut4d_input.ss

#neut4d_input.ss is now an SS-format alignment of all 4d sites in the original alignment

#Step 3. phyloFit
#want to be sure the neutral models are reliable, so run with --init-random, start 5 independent runs of each model (15 total)
#code to run random iterations:

for ITER in 1 2 3 4 5
do
	phyloFit --tree ratiteTree.ver1.nh --init-random --subst-mod SSREV --out-root neut_ver1_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver1_${ITER}.log neut4d_input.ss &> phyloFit_ver1_${ITER}.out &
	phyloFit --tree ratiteTree.ver2.nh --init-random --subst-mod SSREV --out-root neut_ver2_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver2_${ITER}.log neut4d_input.ss &> phyloFit_ver2_${ITER}.out &
	phyloFit --tree ratiteTree.ver3.nh --init-random --subst-mod SSREV --out-root neut_ver3_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver3_${ITER}.log neut4d_input.ss &> phyloFit_ver3_${ITER}.out &
done

#finally, improve all random models to guarantee convergence
for MOD in $(ls neut*.mod);
do
	NEWMOD=${MOD%%.*}
	phyloFit --init-model $MOD --out-root ${NEWMOD}_update --msa-format SS --sym-freqs --log ${MOD}_update.log neut4d_input.ss &> ${MOD}_update.out &
done

#verify convergence of models
#copy each version to _final.mod

cp neut_ver1_1_update.mod neut_ver1_final.mod
cp neut_ver2_2_update.mod neut_ver2_final.mod
cp neut_ver3_2_update.mod neut_ver3_final.mod

#Step 4. Adjust background GC content to be reflective of the average GC content across the non-ancestral genomes in the alignment
###code added to get GC content for each genome, by sampling every 100 bp
mkdir -p baseComp
cd baseComp
for TARGET in $(halStats --genomes /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal)
do
	#output is fraction_of_As fraction_of_Gs fraction_of_Cs fraction_of_Ts
	halStats --baseComp $TARGET,100 /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal > $TARGET.basecomp
done
cd ..

#get average gc content in non-ancestral genomes and update models:
GC=$(cat /n/regal/edwards_lab/ratites/phast/baseComp/??????.basecomp | awk '{SUM+=$2;SUM+=$3;print SUM/42}' | tail -n 1)
for VER in 1 2 3;
do
	modFreqs neut_ver${VER}_final.mod $GC > neut_ver${VER}_corrected.mod
done

#name all ancestral nodes in the tree model
for MOD in 1 2 3
do
	tree_doctor --name-ancestors neut_ver${MOD}_corrected.mod > neut_ver${MOD}_final.named.mod 
done

### END GET NETURAL MODELS ###

### RUN PHYLOP ##
#run halPhyloPMP.py with 12 processors per on each neutral version
#use the _corrected version of each model
for VER in 1 2 3;
do
	mkdir neut_ver$VER
	cp neut_ver${VER}_corrected.mod neut_ver$VER
	cd neut_ver$VER
	halPhyloPMP.py --numProc 12 /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal galGal neut_ver${VER}_corrected.mod galGal_phyloP_ver$VER.wig &> halPhyloP_galGal.log &
	cd ..
done

#also run with ostrich reference
for VER in 1 2 3;
do
	cd neut_ver$VER
	halPhyloPMP.py --numProc 12 /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal strCam neut_ver${VER}_corrected.mod strCam_phyloP_ver$VER.wig &> halPhyloP_strCam.log &
	cd ..
done

#finally also run tree version
for VER in 1 2 3 
do
	cd neut_ver$VER
	halTreePhyloP.py --numProc 24 /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal neut_ver${VER}_corrected.mod . &> halTreePhyloP.log &
	cd ..
done

## END RUN PHYLOP -- NOTE THESE RESULTS NEED TO BE CHECKED ###

### RUNNING PHASTCONS ###
#to run phastCons, we need to take a slightly different approach as there is no direct interface with hal
#so the first step is to export the MAFs that we want, in this case starting with two: chicken, ostrich
#for each MAF, we then run phastCons
#sources: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way and http://compgen.cshl.edu/phast/phastCons-HOWTO.html

#Step 1: Get MAFS
for TARGET in galGal strCam
do
	mkdir -p $TARGET
	cd $TARGET
	hal2mafMP.py --numProc 36 --splitBySequence --sliceSize 5000000 --smallSize 500000 --refGenome $TARGET --noAncestors /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal ${TARGET}_ref.maf &
	cd ..
done

#Step 2. Preprocess MAFS
#fix MAFs
for MAF in $(ls galGal_ref*.maf);
do
	sed -i -e 2d $MAF &
done

#filter duplicates using mafTools
#first step -- remove chicken lines that are from scaffolds other than target with perl script
#second step -- filter with mafDuplicateFilter
for MAF in $(ls galGal_ref_NC*.maf);
do
	./keep_ref_only.pl $MAF &
done

for MAF in $(ls galGal_ref_NC*.temp.maf);
do
	mafDuplicateFilter --maf ${MAF%.temp.maf}.temp.maf > ${MAF%.maf}.pruned.maf &
done

#the output MAFs from mafDuplicateFilter do not guarantee correct strand or order, particularly for galGal specific duplications
#the following code updates / fixes that, first by correcting strand and then order
for MAF in $(ls galGal_ref_NC*.temp.pruned.maf);
do
	mafStrander --maf $MAF --seq galGal --strand + > ${MAF%.pruned.maf}.strand.maf &
done

#finally, we sort the MAF
for FILE in $(ls galGal_ref_NC*.temp.strand.maf);
do
	CHR1=${FILE#galGal_ref_}
	CHR=${CHR1%.temp.strand.maf}
	echo "Processing $CHR"
	mafSorter --maf $FILE --seq galGal.$CHR > $CHR.final.maf &
done

#get rid of temp mafs
rm *.temp*.maf

#split reference into separate files for each chr with samtools
samtools faidx galGal.fa
for FILE in *.final.maf
do
	CHR=${FILE%.final.maf}
	samtools faidx galGal.fa $CHR > $CHR.fa &
done

#Split alignments into chunks
mkdir -p chunks            # put fragments here
for FILE in *.final.maf
do
	CHR=${FILE%.final.maf}
	msa_split $FILE --in-format MAF --refseq $CHR.fa --windows 1000000,0 --out-root chunks/$CHR --out-format SS --min-informative 1000 --between-blocks 5000 2> $CHR.split.log &
done

#Step 3. Estimate rho for chunks and do preliminary phastCons runs
#This will be done with slurm, processing batches of 10 alignments each

#set up job array input
mkdir -p rho
cd rho
ls ../chunks > files
split -a 3 -d -l 10 files part. #make file parts
sbatch est_rho.sh

#run local rerun
./est_rho_local.sh

#kill processes that have not converged after 2 days
#Next -- average rho to get a global rho estimate

ls mods/*.cons.mod > cons.txt
phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
ls mods/*.noncons.mod > noncons.txt
phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod 

#Next -- run phastCons to predict conserved elements on each target segment
sbatch run_phastCons.sh

#Next -- merge predictions and estimate coverage, look at length, other tuning measures
#merge predictions
cat ELEMENTS/*.bed > initCNEEs.bed
cp /n/regal/edwards_lab/phil/PseudoSearch/Final/Chicken_CDS_HLO.bed  chicken_genes.bed
cp /n/regal/edwards_lab/tsackton/ratites/phast/accel/LoweCNEEs.galGal4.bed.fixed lowe_cnees.bed

bedtools coverage -a chicken_genes.bed -b initCNEEs.bed > gene.coverage
awk '{len += $9 ; cov += $8} END {print cov/len}' gene.coverage

#ite1 = 0.558639

bedtools coverage -a lowe_cnees.bed -b initCNEEs.bed > lowe.coverage
awk '{len += $7 ; cov += $6} END {print cov/len}' lowe.coverage

#ite1 = 0.725789

#these are a bit low, so going to try again
#it also seems like perhaps the cnees are a bit fragmentary

#also get total length of elements
awk '{sum+=$3-$2}END{print sum}' phastCons/ite1/initCNEEs.bed 

#83094900

grep "^NC" galGal.chromsizes | awk '{sum+=$2}END{print sum}' 

#1004818361

#implies ~8.2% conserved -- slightly higher than in the feather paper, but probably not surprising given the bird bias in our alignment

#Step 4. Parameter tuning and sensitivity analysis
#going to do this by sampling a subset of mafs to calculate rho, averaging the models, and computing with the averaged models

#ite1 is targetcoverage = 0.3, length = 45
#try target coverage from 0.2 to 0.4, length from 20 to 70

cd ite2
./iterate_phastCons.sh 0.40 20
cd ..
cd ite3
./iterate_phastCons.sh 0.40 45
cd ..
cd ite4
./iterate_phastCons.sh 0.40 70
cd ..
cd ite5
./iterate_phastCons.sh 0.20 20
cd ..
cd ite6
./iterate_phastCons.sh 0.20 45
cd ..
cd ite7
./iterate_phastCons.sh 0.20 70
cd ..
cd ite8
./iterate_phastCons.sh 0.30 70
cd ..
cd ite9
./iterate_phastCons.sh 0.30 20

#get predictions for each iteration
for i in 1 2 3 4 5 6 7 8 9
do
	echo "Processing iteration $i"
	cat ite$i/ELEMENTS/*.bed > cons_ite${i}.bed
	bedtools coverage -a chicken_genes.bed -b cons_ite${i}.bed > gene.coverage.ite${i}
	bedtools coverage -a lowe_cnees.bed -b cons_ite${i}.bed  > lowe.coverage.ite${i}		
	awk '{len += $9 ; cov += $8} END {print cov/len}' gene.coverage.ite${i}
	awk '{len += $7 ; cov += $6} END {print cov/len}' lowe.coverage.ite${i}
	awk '{sum+=$3-$2}END{print sum}' cons_ite${i}.bed
done

#get counts of elements per exon
for i in 1 2 3 4 5 6 7 8 9 
do 
	bedtools intersect -a chicken_genes.bed -b cons_ite${i}.bed -c > ite${i}.gene.count
	bedtools intersect -a lowe_cnees.bed -b cons_ite${i}.bed -c > ite${i}.lowe_cnee.count
done

#Step 5. FINAL PHASTCONS RUNS#

#First need to reprocess the MAFs
#going back to get MAFs for all the small bits -- same code as above but with NW instead of NC
#second step -- filter with mafDuplicateFilter
for MAF in $(ls galGal.NW*.maf);
do
	../keep_ref_only.pl $MAF &
done

for MAF in $(ls galGal_ref_NW*.temp.maf);
do
	mafDuplicateFilter --maf ${MAF%.temp.maf}.temp.maf > ${MAF%.maf}.pruned.maf &
done

#the output MAFs from mafDuplicateFilter do not guarantee correct strand or order, particularly for galGal specific duplications
#the following code updates / fixes that, first by correcting strand and then order
for MAF in $(ls galGal_ref_NW*.temp.pruned.maf);
do
	mafStrander --maf $MAF --seq galGal --strand + > ${MAF%.pruned.maf}.strand.maf &
done

#finally, we sort the MAF
for FILE in $(ls galGal_ref_NW*.temp.strand.maf);
do
	CHR1=${FILE#galGal_ref_}
	CHR=${CHR1%.temp.strand.maf}
	echo "Processing $CHR"
	mafSorter --maf $FILE --seq galGal.$CHR > $CHR.final.maf &
done

#get rid of temp mafs
rm *.temp*.maf

#all mafs are now in final_mafs directory
#need to split the larger chromosomes since phastCons seems to die on chromosomes bigger than ~50 MB
#going to use overlapping 40 MB windows to guarantee no conserved elements fall in breakpoints of windows
#get chrs again, this time including the small bits
#split reference into separate files for each chr with samtools
cd final_mafs
cp ../galGal.fa .
cp ../galGal.fa.fai .
for FILE in *.final.maf
do
	CHR=${FILE%.final.maf}
	samtools faidx galGal.fa $CHR > $CHR.fa
done
rm galGal.fa
rm galGal.fa.fai

#estimate rho using same chunks as before, but this time on the whole data set
sbatch est_rho.sh 1
sbatch est_rho.sh 2

#give runtime of 2 days -- anything not finished will be killed and ignored in the averaging below
 
for VER in 1 2
do
	cd rho_${VER}
	ls mods/*.cons.mod > cons.txt
	phyloBoot --read-mods '*cons.txt' --output-average ave.cons.mod 
	ls mods/*.noncons.mod > noncons.txt
	phyloBoot --read-mods '*noncons.txt' --output-average ave.noncons.mod 
	cd ..
done

#MAKE FINAL SS ALIGNMENTS FOR PHASTCONS

#split and make ratite-specific alignments
mkdir -p split_all
mkdir -p split_nr
for FILE in *.final.maf
do
	CHR=${FILE%.final.maf}
	./make_inputs.sh $CHR $FILE &
done

##make_inputs.sh file:
#!/bin/bash
FILE=$2
CHR=$1
#make chr SS
msa_view $FILE --out-format SS --refseq $CHR.fa --unmask 1> $CHR.ss 2> $CHR.convert.log
#split chr SS
msa_split $CHR.ss --in-format SS --refseq $CHR.fa --windows 40000000,5000000 --out-root ./split_all/$CHR --out-format SS --between-blocks 500000 2> $CHR.split.log
#extract ratite-only alignment
for ALLF in $(ls split_all/$CHR*)
do
	SAMP=${ALLF##*/}
	msa_view $ALLF --in-format SS --seqs rhePen,rheAme,strCam,aptHaa,aptRow,aptOwe,casCas,droNov --exclude --out-format SS 1> ./split_nr/$SAMP 2> $CHR.remove.log
done

#finally, run phastCons with the estimated rho models on the final mafs
./run_phastCons_local.sh 1
./run_phastCons_local.sh 2
./run_phastCons_noRatite.sh 2

##PROCESSING FINAL CONSERVED ELEMENT BEDS###

#merge bed files and remove duplicate lines sort -k1,1 -k2,2n -k3,3n -k6,6 -u 
cat final_run_ver2/ELEMENTS/*.bed | perl -pe 's/^(\w+\.\d)\.\S+/$1/' | sort -k1,1 -k2,2n -k3,3n -k6,6 -u > most_conserved_tree2.bed
#verify no overlapping intervals
bedtools merge -i most_conserved_tree2.bed -c 4 -o collapse -delim "," | grep -c ","

cat final_run_ver1/ELEMENTS/*.bed | perl -pe 's/^(\w+\.\d)\.\S+/$1/' | sort -k1,1 -k2,2n -k3,3n -k6,6 -u > most_conserved_tree1.bed
bedtools merge -i most_conserved_tree1.bed -c 4 -o collapse -delim "," | grep -c "," 

#need to fix the chromosomes in the no-ratite version because of the source files
cat final_run_ver2_noratite/ELEMENTS/*.bed | perl -pe 's/^(\w+\.\d)\.\S+/$1/' | sort -k1,1 -k2,2n -k3,3n -k6,6 -u > most_conserved_noratite_tree2.bed
bedtools merge -i most_conserved_noratite_tree2.bed -c 4 -o collapse -delim "," | grep -c ","

#rename / renumber
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="ce"NR; print}' most_conserved_tree1.bed > most_conserved_final.tree1.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="ce"NR; print}' most_conserved_tree2.bed > most_conserved_final.tree2.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {$4="nrce"NR; print}' most_conserved_noratite_tree2.bed > most_conserved_noratite_final.tree2.bed

##DONE##

##most_conserved_final.tree1.bed and most_conserved_final.tree2.bed contain the final bed elements
##now we'll post-process them with new code