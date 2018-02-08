#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

#pass the list of loci
my $list= shift(@ARGV);
chomp($list);

#store the batch & loci to run...
my %files;
open(LOCI, "<$list") || die "Cannot open $list for reading: $!\n";
while(<LOCI>){
	chomp($_);
	my @spla= split(/\t/, $_);
	#store keyed by locus name, with batch as value
	$files{$spla[1]}= $spla[0];
}
my $scal= scalar(keys %files);
print "Stored list of: $scal loci to align with prank\n";

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart)
my $baselist= basename($list);
my $logfile= 'prank_logfile_' . $baselist;
my %done;
if (-e $logfile){
	open(DONEJOBS, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONEJOBS>){
		if ($_=~ /^(.+?)(\s)(done)$/){
			$done{$1}= 'exists';
		}
	}
	close(DONEJOBS);

	#reopen for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

#create a top-level directory to hold output, unless it already exists
my $topdir= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/PAML_aligned_PRANK_fastas/';
unless (-e $topdir){
	mkdir $topdir || die "Cannot create $topdir: $!\n";
}

for my $locus (sort {$a cmp $b} (keys %files)){
	print "Processing: $locus\n";
	my $batch= $files{$locus};

	#if we've already aligned this locus...skip
	if (exists $done{$locus}){
		next;
	}
	#create a subdirectory for this batch...
	my $outdir= $topdir . $batch . '/';
	unless (-e $outdir){
		mkdir $outdir || die "Cannot create $outdir: $!\n";
	}
	#full path to input file
	my $infile= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/PAML_input_PRANK_fastas/' . $batch . '/' . $locus;
	#full path to output file
	my $outfile= $outdir . $locus;

	#command for prank alignment...
	my $cmd= "prank -d=$infile -o=$outfile -codon";
	
	#run prank as system command within perl wrapper...
	system($cmd);
	
	#record completion in logfile
	print LOG "$locus done\n";
}
