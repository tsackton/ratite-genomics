#phast code for runs with moa

#starting point is 1) 4d alignment with moa added, and 2) cnee alignment with moa added
#goal is to reestimate neutral model and rerun phyloP 

#neutral model: tree files are ratiteTree_withMoa.ver?.nh

#get species
SPECIES=$(nw_labels ratiteTree_withMoa.ver1.nh | sort | tr '\n' ',')
#get alignment
MSAFILES=allspecies_4d_concatenated.fasta
#convert to SS

#Note: this includes soft-masked sites since they were converted to uppercase when extracting from MAF

msa_view --aggregate ${SPECIES%?} --in-format FASTA --out-format SS --unordered-ss $MSAFILES > neut4d_input_withMoa.ss

#neut4d_input_withMoa.ss is now an SS-format alignment of all 4d sites in the original alignment

#now get neutral models
#want to be sure the neutral models are reliable, so run with --init-random, start 5 independent runs of each model (15 total)
#code to run random iterations:

for ITER in 1 2 3 4 5
do
	phyloFit --tree ratiteTree_withMoa.ver1.nh --init-random --subst-mod SSREV --out-root neut_ver1_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver1_${ITER}.log neut4d_input_withMoa.ss &> phyloFit_ver1_${ITER}.out &
	phyloFit --tree ratiteTree_withMoa.ver2.nh --init-random --subst-mod SSREV --out-root neut_ver2_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver2_${ITER}.log neut4d_input_withMoa.ss &> phyloFit_ver2_${ITER}.out &
	phyloFit --tree ratiteTree_withMoa.ver3.nh --init-random --subst-mod SSREV --out-root neut_ver3_${ITER} --msa-format SS --sym-freqs --log phyloFit_ver3_${ITER}.log neut4d_input_withMoa.ss &> phyloFit_ver3_${ITER}.out &
done

#finally, improve all random models to guarantee convergence
for MOD in $(ls neut*.mod);
do
	NEWMOD=${MOD%%.*}
	phyloFit --init-model $MOD --out-root ${NEWMOD}_update --msa-format SS --sym-freqs --log ${MOD}_update.log neut4d_input_withMoa.ss &> ${MOD}_update.out &
done

#final versions of neutral models
cp neut_ver1_4_update.mod neut_ver1_final.mod
cp neut_ver2_4_update.mod neut_ver2_final.mod
cp neut_ver3_4_update.mod neut_ver3_final.mod

#get average gc content in non-ancestral genomes and update models:
GC=$(cat /n/regal/edwards_lab/ratites/wga/phast/baseComp/??????.basecomp | awk '{SUM+=$2;SUM+=$3;print SUM/42}' | tail -n 1)
for VER in 1 2 3;
do
	modFreqs neut_ver${VER}_final.mod $GC > neut_ver${VER}_corrected.mod
done

#name all ancestral nodes in the tree model
for MOD in 1 2 3
do
	tree_doctor --name-ancestors neut_ver${MOD}_corrected.mod > neut_ver${MOD}_final.named.mod 
done

#get rho from previous run of phastCons
cp ../../phastCons/ave.*.mod .
#rho = 0.3168

#fix partition file to be a bedfile
awk 'BEGIN{FS="[ -]"; OFS="\t"}{print "ref", $4-1, $5, $2}' allspecies_cnee_concat_partitions > cnees.bed

