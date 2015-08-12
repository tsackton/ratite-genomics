#SNAP training (droNov single):

cd $HOME/ratite_scratch/training/droNov

#Extract good gene models from maker GFF
maker2zff ~/ratite_scratch/maker/dronov/droNov.genome.orig.gff #using default parameters

#check quality of zff files
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1

#split into training set 1, training set 2, and testing set
fathom -split -number 3 uni.ann uni.dna
mv 0.ann training1.ann
mv 0.dna training1.dna
mv 1.ann training2.ann
mv 1.dna training2.dna
cat training?.ann > training_droNov.ann
cat training?.dna > training_droNov.dna
mv 2.ann testing.ann
mv 2.dna testing.dna
fathom training_droNov.ann training_droNov.dna -export 1000 -plus > uni-plus.log 2>&1

#train
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl emu params/ > emu.hmm

#SNAP training (notPer single):

cd $HOME/ratite_scratch/training/notPer

#Extract good gene models from maker GFF
maker2zff ~/ratite_scratch/maker/notper/notPer.genome.orig.gff #using default parameters

#check quality of zff files
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1

#split into training set 1, training set 2, and testing set
fathom -split -number 3 uni.ann uni.dna
mv 0.ann training1.ann
mv 0.dna training1.dna
mv 1.ann training2.ann
mv 1.dna training2.dna
cat training?.ann > training_notPer.ann
cat training?.dna > training_notPer.dna
mv 2.ann testing.ann
mv 2.dna testing.dna
fathom training_notPer.ann training_notPer.dna -export 1000 -plus > uni-plus.log 2>&1

#train
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl tinamou params/ > tinamou.hmm

#SNAP training (strCam single):

cd $HOME/ratite_scratch/training/strCam

#get data from NCBI
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Struthio_camelus_australis/GFF/ref_ASM69896v1_top_level.gff3.gz
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000698965.1_ASM69896v1/GCF_000698965.1_ASM69896v1_genomic.fna.gz
gunzip *.gz

#filter GFF to keep just best models
module load cufflinks
gffread ref_ASM69896v1_top_level.gff3 -g GCF_000698965.1_ASM69896v1_genomic.fna -C -J --no-pseudo -o strCam.good.gff
echo -e "##FASTA" >> strCam.good.gff 
cat GCF_000698965.1_ASM69896v1_genomic.fna >> strCam.good.gff

#maker2zff to parse filtered GFF
maker2zff -n strCam.good.gff

#check quality of zff files
fathom genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
fathom genome.ann genome.dna -validate > validate.log 2>&1
fathom genome.ann genome.dna -categorize 1000 > categorize.log 2>&1

#split into training set 1, training set 2, and testing set
fathom -split -number 3 uni.ann uni.dna
mv 0.ann training1.ann
mv 0.dna training1.dna
mv 1.ann training2.ann
mv 1.dna training2.dna
cat training?.ann > training_strCam.ann
cat training?.dna > training_strCam.dna
mv 2.ann testing.ann
mv 2.dna testing.dna
fathom training_strCam.ann training_strCam.dna -export 1000 -plus > uni-plus.log 2>&1

#train
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl ostrich params/ > ostrich.hmm

#SNAP training (ratite merged):

cd $HOME/ratite_scratch/training/ratite
cp ../strCam/training1.ann strcam.ann
cp ../strCam/training1.dna strcam.dna
cp ../droNov/training1.ann dronov.ann
cp ../droNov/training1.dna dronov.dna

cat strcam.ann dronov.ann > training_ratite.ann
cat strcam.dna dronov.dna > training_ratite.dna
fathom training_ratite.ann training_ratite.dna -export 1000 -plus > uni-plus.log 2>&1

#train
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
hmm-assembler.pl ratite params/ > ratite.hmm

