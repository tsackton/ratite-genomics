module load cufflinks

for SP in aptHaa aptRow aptOwe casCas droNov cryCin eudEle notPer rheAme rhePen
do
	mkdir -p $SP
	gffread --no-pseudo -C -T -o $SP.gtf $SP.genome.gff
	grep -P "\tCDS\t" $SP.gtf > ${SP}_filt.gtf
	gtfToGenePred ${SP}_filt.gtf ${SP}.gp
	genePredToBed ${SP}.gp ${SP}.bed
	mv $SP.bed $SP
	rm $SP.*
done

for SP in $(ls alignability/); do mkdir -p phyloP/$SP; twoBitInfo twoBit/$SP.2bit chrom.size; wigToBigWig ../../phast/phyloP/neut_ver2/${SP}_phyloP.wig chrom.size phyloP/$SP/$SP.bw; done

for SP in aptHaa aptRow aptOwe casCas droNov cryCin eudEle notPer rheAme rhePen
do
	cd $SP
	twoBitInfo ../../twoBit/$SP.2bit chrom.size
	sort -k1,1 -k2,2n ~/ratite_store/annotation/maker_final/$SP.bed > input.bed
	bedToBigBed -type=bed12 -extraIndex=name input.bed chrom.size $SP.bb
	rm chrom.size
	rm input.bed
	cd ..
done
	
hal2assemblyHub.py /n/regal/edwards_lab/ratites/wga/ratite_final_20150627/ratiteAlign.hal ratiteHub \
	--hub=ratiteHub \
	--shortLabel="RatiteHub" \
	--longLabel="Comparative assembly hub for ratite genomes" \
	--email=tsackton@oeb.harvard.edu \
	--rename=common_names.txt \
	--lodDir=/n/regal/edwards_lab/ratites/assemblyHub/inputs/lodFinal \
	--lodTxtFile=/n/regal/edwards_lab/ratites/assemblyHub/inputs/lod.txt \	
	--twobitdir=/n/regal/edwards_lab/ratites/assemblyHub/inputs/twoBit \
	--bedDirs=/n/regal/edwards_lab/ratites/assemblyHub/inputs/Genes \
	--finalBigwigDirs=/n/regal/edwards_lab/ratites/assemblyHub/inputs/GC,/n/regal/edwards_lab/ratites/assemblyHub/inputs/alignability,/n/regal/edwards_lab/ratites/assemblyHub/inputs/phyloP \
	--finalBigBedDirs=/n/regal/edwards_lab/ratites/assemblyHub/inputs/MAKER \
	--maxThreads 24 \
