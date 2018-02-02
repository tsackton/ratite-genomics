#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

#pass the list of loci to align (e.g. can split cnee_alignment_batches to run in parallel)
my $list= shift(@ARGV);
my %toalign;
open(LIST, "<$list") || die "Cannot open $list for reading: $!\n";
while(<LIST>){
	chomp($_);
	my @split= split(/\t/, $_);
	#store, keyed by locus name, & with batch as value
	$toalign{$split[1]}= $split[0];
}
close(LIST);
my $scalar= scalar(keys %toalign);
print "\nPreparing to align: $scalar loci\n";

my $outdir= '/n/regal/edwards_lab/acloutier/moa_cnees/input_fastas_allspecies_cnees_aligned/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}
my $processed= 0;
for my $locus (sort {$a cmp $b} (keys %toalign)){
	my $batch= $toalign{$locus};

	my $subdir= $outdir . $batch . '/';
	unless (-e $subdir){
		mkdir $subdir || die "Cannot create $subdir: $!\n";
	}
	my $fasta_in= '/n/regal/edwards_lab/acloutier/moa_cnees/input_fastas_allspecies_cnees/' . $batch . '/' . $locus . '.fasta';
	my $fasta_out= $subdir . $locus . '_aln.fasta';
	#skip if already done (use file test to check for output file of nonzero size)
	if (-s $fasta_out){
		#print "Skipping: $locus (already done)\n";
		next;
	}
		
	#otherwise, run/re-run mafft...
	#print "Run/rerun mafft on: $locus\n";
	
	#command to run mafft
	my $cmd= 'ginsi ' . $fasta_in . ' > ' . $fasta_out;
	system($cmd);
	$processed++;
	if ($processed=~ /^(\d+)(0{2})$/){
		print "Processed: $processed\n";
	}
}
