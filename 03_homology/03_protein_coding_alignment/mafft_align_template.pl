#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

#pass the batch number...
my $batch= shift(@ARGV);
chomp($batch);

#glob all of the unaligned fastas for input batch
my @glob= </n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas/$batch/*.fa>;

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart on Odyssey cluster)
my $logfile= 'mafft_align_logfile_' . $batch;
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

my $outdir= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned/' . $batch . '/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}

for my $filepath (sort {$a cmp $b} (@glob)){
	my $basename= basename($filepath);
	my $locus;
	if ($basename=~ /^(.+)(\.fa)$/){
		$locus= $1;
	}
	else{
		print "Error parsing locus from: $basename\n";
		die;
	}

	#if we've already aligned this locus...skip
	if (exists $done{$locus}){
		next;
	}
	
	my $fasta_out= $outdir . $locus . '_aln.fa';

	#command to run mafft
	my $cmd= 'linsi ' . $filepath . ' > ' . $fasta_out;
	system($cmd);

	print LOG "$locus done\n";	
}
