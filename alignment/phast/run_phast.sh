#make tree
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

#NOTE: 9/22/2015 THIS CODE DOES NOT WORK DUE TO HAL2MAF BUG/UNDOCUMENTED FEATURE WITH REFSEQUENCE / REFTARGETS
#hal2mafMP.py --numProc 48 --refGenome galGal --noAncestors --noDupes --refTargets galGal4_4d.bed /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal neut4d_input_galGal_ref.maf

#NEW CODE BELOW:
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
msa_view --aggregate ${SPECIES%?} --in-format MAF --out-format SS --unordered-ss $MSAFILES > neut4d_input.ss

#NEED TO RERUN EVERYTHING BELOW HERE#


#run phyloFit
#want to be sure the neutral models are reliable, so run with --init-random, start 5 independent runs of each model (15 total)
#also because of convergence error on ver1, rerun with --init-model to retry 
#code to run random iterations:

for ITER in 1 2 3 4 5
do
	phyloFit --tree ratiteTree.ver1.nh --init-random --subst-mod SSREV --out-root neut_ver1_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver1_${ITER}.log neut4d_input.ss &> phyloFit_ver1_${ITER}.out &
	phyloFit --tree ratiteTree.ver2.nh --init-random --subst-mod SSREV --out-root neut_ver2_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver2_${ITER}.log neut4d_input.ss &> phyloFit_ver2_${ITER}.out &
	phyloFit --tree ratiteTree.ver3.nh --init-random --subst-mod SSREV --out-root neut_ver3_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver3_${ITER}.log neut4d_input.ss &> phyloFit_ver3_${ITER}.out &
done

#if any models fail to converge as indicated by the log file, rerun with --init-model

#finally, improve all random models with same approach as above
for MOD in $(ls neut*.mod);
do
	NEWMOD=${MOD%%.*}
	phyloFit --init-model $MOD --out-root ${NEWMOD}_update --msa-format SS --sym-freqs --log ${MOD}_update.log neut4d_input.ss &> ${MOD}_update.out &
done

#verify convergence of models
#copy each version to _final.mod
#run halPhyloPMP.py with 12 processors per on each neutral version
for VER in 1 2 3;
do
	mkdir neut_ver$VER
	cp neut_ver${VER}_final.mod neut_ver$VER
	cd neut_ver$VER
	halPhyloPMP.py /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal galGal neut_ver${VER}_final.mod galGal_phyloP_ver$VER.wig &> halPhyloP.log &
	cd ..
done

#also run with ostrich reference
for VER in 1 2 3;
do
	cd neut_ver$VER
	halPhyloPMP.py /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal strCam neut_ver${VER}_final.mod strCam_phyloP_ver$VER.wig &> halPhyloP_strCam.log &
	cd ..
done


	
