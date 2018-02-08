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
my $aln= 'UCEs_mafft_concat.phylip';

#partition file
my $partitions= 'UCEs_mafft_concat_partitions';

#substitution model
my $model= 'GTRGAMMA';

#create a binary input file for ExaML:
print "\nGenerating binary input file...\n\n";
my $model2= 'DNA';
my $examl= "parse-examl -s $aln -q $partitions -m $model2 -n UCEs_mafft_concat";
system($examl);

#generate 20 complete random starting trees:
for (my $i= 1; $i<= 20; $i++){
	print "\n\nGenerating files for run$i\n";
	#random number
	my $range= 100000;
	my $min= 1;
	my $rand1= int(rand($range)) + $min;
	#output file name
	my $name= 'StartingTree' . $i;
	if (-s $name){
		next;
	}
	else{
		#raxml command
		my $raxml= "raxmlHPC-SSE3 -y -d -m $model -s $aln -q $partitions -p $rand1 -n $name";
		system($raxml);
	
		#rename output file to just StartingTree# (so both random & parsimony starting trees can be run with a single slurm array command)
		my $oldfile= "RAxML_randomTree.StartingTree$i";
		my $newfile= "StartingTree$i";
		move($oldfile,$newfile) || die "File renaming failed for $oldfile -> $newfile\n";
	}
}
#plus, do one run with random stepwise addition order parsimony starting tree
for (my $j= 21; $j<= 21; $j++){
	print "\n\nGenerating files for run$j\n";
	#random number
	my $range2= 100000;
	my $min2= 1;
	my $rand2= int(rand($range2)) + $min2;
	#output file name
	my $name2= 'StartingTree' . $j;
	if (-s $name2){
		next;
	}
	else{
		#raxml command
		my $raxml2= "raxmlHPC-SSE3 -y -m $model -s $aln -q $partitions -p $rand2 -n $name2";
		system($raxml2);

		my $oldfile2= "RAxML_parsimonyTree.StartingTree$j";
		my $newfile2= "StartingTree$j";
		move($oldfile2,$newfile2) || die "File renaming failed for $oldfile2 -> $newfile2\n";
	}
}
