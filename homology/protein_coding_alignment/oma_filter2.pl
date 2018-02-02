#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

#pass the list of loci (e.g. batch1_list1)...
#NB- lists can be created using bash commands ls & split...
my $batchlist= shift(@ARGV);
chomp($batchlist);
my $batch;
if ($batchlist=~ /^(batch)(\d+)(_)(list)(\d+)$/){
	$batch= $1 . $2;
}
else{
	print "Error parsing batch from: $batchlist\n";
	die;
}

#store all of the loci from the passed list
open(LIST, "<$batchlist") || die "Cannot open $batch for reading: $!\n";
my %list;
while(<LIST>){
	chomp($_);
	$list{$_}= 'exists';
}
close(LIST);

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart)
my $logfile= 'jarvis_filter_logfile_' . $batchlist;
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

my $jfile_in= 'spotProblematicSeqsBase-W12S4.py';

for my $key1 (sort {$a cmp $b} (keys %list)){
	my $filepath= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned_filter1/' . $batch . '/' . $key1;
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
	print "Processing: $locus\n";
	#create an input file that lists just current locus...
	my $locfile= $locus;
	open(LIST, ">$locfile") || die "Cannot open $locfile 'listfile' for writing: $!\n";
	print LIST "$basename\n";
	close(LIST);	

	#create a new jarvis filtering python script for this input fasta
	my $jfile_out= 'jarvis_filter_' . $basename . '.py';
	open(OUT, ">$jfile_out") || die "Cannot open $jfile_out for writing: $!\n";
	
	#go through the input Jarvis file...
	open(IN, "<$jfile_in") || die "Cannot open $jfile_in for reading: $!\n";
	while(<IN>){
		chomp($_);
		#skip the commented-out fileLocation lines
		if ($_=~ /(\#)(\s*)(fileLocation=)/){
			next;
		}
		#if we've reached the filename rstrip line...
		elsif ($_=~ /(filename=filename)(\.)(rstrip)/){
			#print out line
			print OUT "$_\n";
			#then, print out the file location for current list
			my $fileloc= "    fileLocation=\"/n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned_filter1/$batch/\%s\"\%(filename)";
			print OUT "$fileloc\n";
		}
		#otherwise, print out line...
		else{
			print OUT "$_\n";
		}
	}
	close(IN);
	close(OUT);

	#run jarvis filter as a system command
	my $cmd= "python $jfile_out $locfile";
	system($cmd);

	unlink($jfile_out);
	unlink($locfile);
	#and, delete the outfiles/tempfiles that aren't wanted...
	my $del1= 'W12S1G6_' . $locus . '_WrongC.txt';
	unlink($del1);
	my $del2= 'W12S1G6_' . $locus . '_clean.temp';
	unlink($del2);
	
	print LOG "$locus done\n";
}
