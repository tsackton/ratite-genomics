##Training Augustus## (based on: http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html)
#make GFFs

##NOTE: ../droNov/training_droNov.ann and similar are the files produced for SNAP training##
##path needs to be set properly for this to work##

zff2gff3.pl ../droNov/training_droNov.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > droNov.gff3
zff2gff3.pl ../ratite/training_ratite.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > ratite.gff3
zff2gff3.pl ../notPer/training_notPer.ann | perl -plne 's/\t(\S+)$/\t\.\t$1/' > notPer.gff3

#setup
new_species.pl --species=emu
new_species.pl --species=ratite
new_species.pl --species=tinamou

#make genbank
gff2gbSmallDNA.pl droNov.gff3 ../droNov/training_droNov.dna 1000 droNov.gb
gff2gbSmallDNA.pl notPer.gff3 ../notPer/training_notPer.dna 1000 notPer.gb 
gff2gbSmallDNA.pl ratite.gff3 ../ratite/training_ratite.dna 1000 ratite.gb

#initial etraining
etraining --species=emu droNov.gb
etraining --species=tinamou notPer.gb
etraining --species=ratite ratite.gb

#optimize
optimize_augustus.pl --species=emu droNov.gb --cpus=12 1>emu.stdout 2>emu.stderr &
optimize_augustus.pl --species=ratite ratite.gb --cpus=12 1>ratite.stdout 2>ratite.stderr &
optimize_augustus.pl --species=tinamou notPer.gb --cpus=12 1>tinamou.stdout 2>tinamou.stderr &


