#!/usr/bin/perl
use strict;
use warnings;

use POSIX;

#modules needed on Odyssey cluster
#source new-modules.sh
#module load perl
#module load perl-modules
#module load gcc/4.8.2-fasrc01
#module load raxml/8.1.5-fasrc03
#module load openmpi/1.10.0-fasrc01
#module load ExaML/2.0.4-fasrc02

#alignment file
my $aln= 'CNEEs_mafft_concat.phylip';

#partition file
my $partitions= 'CNEEs_mafft_concat_partitions';

#substitution model
my $model= 'GTRGAMMA';

#clean up any partially created files from serial requeue stop/restart...
my @glob1= <RAxML_randomTree*>;
for my $path1 (@glob1){
	unlink($path1);
}
my @glob2= <RAxML_info.run1_BS*>;
for my $path2 (@glob2){
	unlink($path2);
}
my @glob3= <CNEEs_mafft_concat.phylip.BS*>;
for my $path3 (@glob3){
	unlink($path3);
}
my @glob4= <run1_BS*_aln.binary>;
for my $path4 (@glob4){
	unlink($path4);
}

#create 100 bootstrap datasets
my $range= 100000;
my $min= 1;
my $rand2= int(rand($range)) + $min;
my $raxml_bs= "raxmlHPC-SSE3 -s $aln -q $partitions -f j -N 100 -b $rand2 -m GTRGAMMA -n run1_BS";
system($raxml_bs);

for (my $i= 0; $i<= 199; $i++){
	print "\n\nGenerating files for run1 BS$i\n";
	#random number
	my $rand1= int(rand($range)) + $min;
	#output file name
	my $name= 'run1_BS' . $i;

	#raxml command to generate a random starting tree for current bootstrap replicate
	my $raxml= "raxmlHPC-SSE3 -y -d -m $model -s $aln -q $partitions -p $rand1 -n $name";
	system($raxml);

	#convert bootstrap replicate for this run to ExaML binary format
	my $parseexaml= "parse-examl -s CNEEs_mafft_concat.phylip.BS$i -q $partitions -m DNA -n run1_BS$i\_aln";
	system($parseexaml);
}
