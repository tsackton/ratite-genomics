#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use List::Util qw(shuffle);

#specify run # (for independent full run of mp-est)
my $run= 1;

#create a directory to hold output
my $topdir= 'introns_mafft_mpest/';

unless (-e $topdir){
	mkdir($topdir) || die "Cannot create $topdir: $!\n";
}

#create mpest input & output subdirectories
my $subdir_in= $topdir . 'mpest_input_run' . $run . '/';
#first, remove any existing files (from serial requeue stop/restart)
if (-e $subdir_in){
	my $remove= "rm -r $subdir_in";
	system($remove);
}
#then, make input subdirectory afresh...
mkdir($subdir_in) || die "Cannot create $subdir_in: $!\n";

my $subdir_out= $topdir . 'mpest_output_run' . $run . '/';
unless (-e $subdir_out){
	mkdir($subdir_out) || die "Cannot create $subdir_out: $!\n";
}

#glob all rooted bootstrap replicate files
my @glob= <introns_mafft_rerooted_bootstrap_trees/*.tre>;

my $scalar_glob= scalar(@glob);
print "\nPreparing to construct mpest input files for 500 bootstrap replicates of $scalar_glob loci...\n\n";

#go through each locus bootstrap file & print each bootstrap of each gene to a randomly numbered 'genetree' mpest input file

my $filecount= 0;
for my $filepath (sort {$a cmp $b} (@glob)){
	$filecount++;
	
	#open current file for reading...
	open(INFILE, "<" . $filepath) || die "Cannot open $filepath for reading: $!\n";

	#create an array of 'bootnums' from 1 to tot # of bs replicates...
	my @bootnums;
	for (my $z= 1; $z <= 500; $z++){
		push(@bootnums, $z);
	}
	#shuffle the list of bootnums to grab a random number between 1-500 each time
	#(this will randomize which bootstrap replicates get assigned to each mpest 'genetree' file
	#so we don't generate the exact same genetree files for multiple independent runs of mpest)
	my @bootnums2= shuffle(@bootnums);
	my $treecount2= 0;
	while(<INFILE>){
		chomp($_);
		$treecount2++;
		my $treecount= shift(@bootnums2);
		my $useline;
		#remove the [&R] from the start of the line...
		if ($_=~ /^(\[&R\]\s)(.+)$/){
			$useline= $2;
		}
		else{
			print "Error parsing: '[&R] ' from start of line...\n";
			die;
		}
		#open output genetree file for writing/appending
		my $treefile= $subdir_in . 'genetree' . $treecount . '.tre';
		if (-e $treefile){
			open(TREE, ">>" . $treefile) || die "Cannot open $treefile for appending: $!\n";
		}
		else{
			open(TREE, ">" . $treefile) || die "Cannot open $treefile for writing: $!\n";
		}
		print TREE "$useline\n";
		close(TREE);
	}
	close(INFILE);
	#make sure that we really had 500 replicates for this locus...
	unless($treecount2 == 500){
		print "Error: $filepath had $treecount2 bootstrap replicates (expect 500)\n";
		die;
	}
	if ($filecount=~ /^(\d+)(0{2})$/){
		print "Processed: $filecount\n";
	}
}
print "Finished storing bootstrap replicates for: $filecount input loci...\n\n";

#list of all taxon IDs
my @taxa= qw( anoDid aptHaa aptMan aptOwe aptRow casCas cryCin droNov eudEle galGal notPer rheAme rhePen strCam tinGut );

#map array entries as hash keys
my %taxa= map { $_ => 1} @taxa;

my $scalar_taxa= scalar(keys %taxa);
print "\nFinished storing taxon IDs for: $scalar_taxa taxa...\n\n";

#create an mpest control file for each 'genetree' (1 'genetree' per bootstrap)
print "\nCreating control files to run MP-EST...\n\n";
for (my $index1= 1; $index1 <= 500; $index1++){	
	my $controlfile= $subdir_in . 'control_' . $index1;
	open(CONTROL, ">" . $controlfile) || die "Cannot open $controlfile for writing: $!\n";
	#print everything needed to the control file...
	print CONTROL "genetree$index1\.tre\n";
	print CONTROL "0\n";
	#generate a random number to use as seed
	my $range= 100000;
	my $minimum= 1;
	my $random_number= int(rand($range)) + $minimum;
	print CONTROL "$random_number\n";
	print CONTROL "10\n";
	print CONTROL "$scalar_glob $scalar_taxa\n";
	for my $taxkey (sort {$a cmp $b} (keys %taxa)){
		print CONTROL "$taxkey\t1\t$taxkey\n";
	}
	print CONTROL "0\n";
	close(CONTROL);
}
