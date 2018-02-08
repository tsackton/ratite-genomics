#!/usr/bin/perl
use strict;
use warnings;

use File::Copy;

#store batch info. for all loci
my $batchlist= 'cnee_alignment_batches';
my %batches;
open(BATCHES, "<$batchlist") || die "Cannot open $batchlist for reading: $!\n";
while(<BATCHES>){
	chomp($_);
	if ($_=~ /^(Batch)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		#store keyed by locus
		$batches{$split[1]}= $split[0];
	}
}
close(BATCHES);

#top-level directory to hold parsed output
#NB- comment/uncomment lines depending on which dataset is being processed
my $outdir= 'cnee_branchlengths_ver1/';
#my $outdir= 'cnee_branchlengths_ver2/';
#my $outdir= 'cnee_branchlengths_ver3/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}

#output file to hold model parameters for all CNEEs
my $outfile= $outdir . 'RAxML_info.cnee_branchlengths_ver1';
#my $outfile= $outdir . 'RAxML_info.cnee_branchlengths_ver2';
#my $outfile= $outdir . 'RAxML_info.cnee_branchlengths_ver3';
open(OUT, ">$outfile") || die "Cannoto open $outfile for writing: $!\n";
print OUT "Batch\tLocus\tAlignment_patterns\talpha\tTree-Length\trateA<->C\trateA<->G\trateA<->T\trateC<->G\trateC<->T\trateG<->T\tfreq_pi(A)\tfreq_pi(C)\tfreq_pi(G)\tfreq_pi(T)\n";

my %partitions;
my $total= 0;
my %check;
#for each RAxML run...
for (my $i= 1; $i <= 2254; $i++){
	print "\nProcessing run$i\n";
	my $count= 0;
	#directory holding RAxML output...
	my $raxdir= 'cnee_branchlengths_ver1_raw/' . $i . '/';
	#my $raxdir= 'cnee_branchlengths_ver2_raw/' . $i . '/';
	#my $raxdir= 'cnee_branchlengths_ver3_raw/' . $i . '/';

	#open info file for parsing
	my $raxinfo= $raxdir . 'RAxML_info.cnee_branchlengths_ver1_' . $i;
	#my $raxinfo= $raxdir . 'RAxML_info.cnee_branchlengths_ver2_' . $i;
	#my $raxinfo= $raxdir . 'RAxML_info.cnee_branchlengths_ver3_' . $i;
	unless(-s $raxinfo){
		print "Check run$i (no RAxML_info)\n";
	}
	open(INFO, "<$raxinfo") || die "Cannot open $raxinfo for reading: $!\n";
	
	my $storedlocus;
	my $storedpatterns;
	while(<INFO>){
		chomp($_);
		#if we've reached a line with info. about # alignment patterns...
		if ($_=~ /^(Alignment\sPatterns:\s)(\d+)$/){
			$storedpatterns= $2;
		}
		#if we've reached a line with a locus name & we have stored patterns...
		elsif (($_=~ /^(Name:\smCE)(.+)$/) && (defined($storedpatterns))){
			$partitions{$2}{'patterns'}= $storedpatterns;
			$storedpatterns= '';
		}
		#if we've reached the start of a section of model parameters...
		elsif ($_=~ /^(Model\sParameters\sof\sPartition\s)(\d+?)(\,)(\sName:\s)(.+?)(\,)(\sType\sof\sData:\sDNA)$/){
			my $part= $2;
			my $temploc= $5;
			#parse just the numerical portion to use as key (for better sorting...)
			if ($temploc=~ /^(mCE)(\d+)$/){
				$storedlocus= $2;
			}
			else{
				print "Error parsing numerical portion of: $temploc\n";
				die;
			}
			$partitions{$storedlocus}{'name'}= $temploc;
			$partitions{$storedlocus}{'partition'}= $part;
			#& also store the 'run'...
			$partitions{$storedlocus}{'run'}= $i;
			
			$count++;
			$total++;
			$check{$temploc}= 'exists';
		}
		#then, for each following parameter...
		elsif (($_=~ /^(alpha:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'alpha'}= $2;
		}
		elsif (($_=~ /^(Tree-Length:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'treelength'}= $2;
		}
		elsif (($_=~ /^(rate\sA\s<->\sC:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'AC'}= $2;
		}
		elsif (($_=~ /^(rate\sA\s<->\sG:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'AG'}= $2;
		}
		elsif (($_=~ /^(rate\sA\s<->\sT:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'AT'}= $2;
		}
		elsif (($_=~ /^(rate\sC\s<->\sG:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'CG'}= $2;
		}
		elsif (($_=~ /^(rate\sC\s<->\sT:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'CT'}= $2;
		}
		elsif (($_=~ /^(rate\sG\s<->\sT:\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'GT'}= $2;
		}
		elsif (($_=~ /^(freq\spi\(A\):\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'freqA'}= $2;
		}
		elsif (($_=~ /^(freq\spi\(C\):\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'freqC'}= $2;
		}
		elsif (($_=~ /^(freq\spi\(G\):\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'freqG'}= $2;
		}
		elsif (($_=~ /^(freq\spi\(T\):\s)(.+)$/) && (defined($storedlocus))){
			$partitions{$storedlocus}{'freqT'}= $2;
			
			#this is the last value, so reset storedlocus
			$storedlocus= '';
		}
		else{
			next;
		}
	}
	close(INFO);
}
print "Finished parsing info. for: $total loci total\n";	
my $scalarloci= scalar(keys %partitions);
print "Stored info. for: $scalarloci loci\n\n";
unless ($scalarloci == 284001){
	print "Error!  Expect 284001 loci!\n";
	die;
}

print "Outputting data...\n";
my $output= 0;
for my $key1 (sort {$a <=> $b} (keys %partitions)){
	#get the full name
	my $locus= $partitions{$key1}{'name'};
	#get the batch
	my $batch= $batches{$locus};
	#print out all values
	print OUT "$batch\t$locus\t$partitions{$key1}{'patterns'}\t$partitions{$key1}{'alpha'}\t$partitions{$key1}{'treelength'}\t$partitions{$key1}{'AC'}\t$partitions{$key1}{'AG'}\t$partitions{$key1}{'AT'}\t$partitions{$key1}{'CG'}\t$partitions{$key1}{'CT'}\t$partitions{$key1}{'GT'}\t$partitions{$key1}{'freqA'}\t$partitions{$key1}{'freqC'}\t$partitions{$key1}{'freqG'}\t$partitions{$key1}{'freqT'}\n";

	#copy all of the individual partition trees with estimated branchlengths...
	my $partnum= $partitions{$key1}{'partition'};
	my $runnum= $partitions{$key1}{'run'};
	my $subdir= $outdir . $batch . '/';
	unless (-e $subdir){
		mkdir $subdir || die "Cannot create $subdir: $!\n";
	}
	my $raxdir2= '/n/regal/edwards_lab/acloutier/moa_cnees/cnee_branchlengths_ver1_raw/' . $runnum . '/';
	my $tree1= $raxdir2 . 'RAxML_result.cnee_branchlengths_ver1_' . $runnum . '.PARTITION.' . $partnum;
	my $tree2= $subdir . 'RAxML_result.' . $locus . '_ver1';
	if (-s $tree1){
		copy($tree1,$tree2) || die "Copy failed for: $locus: $!\n";
	}
	else{
		print "Error!  No treefile: $tree1\n";
	}
	$output++;
	if ($output=~ /^(\d+)(0{3})$/){
		print "Output: $output\n";
	}
}
close(OUT);
