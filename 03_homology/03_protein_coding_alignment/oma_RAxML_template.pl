#!/usr/bin/perl
use strict;
use warnings;

use POSIX;
use File::Spec;
use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;

#pass the list of loci
my $list= shift(@ARGV);
chomp($list);

#store the batch & loci to run...
my %files;
open(LOCI, "<$list") || die "Cannot open $list for reading: $!\n";
while(<LOCI>){
	chomp($_);
	my @spla= split(/\t/, $_);
	#store keyed by batch, then full path to locus
	$files{$spla[0]}{$spla[1]}= 'exists';
}

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart)
my $baselist= basename($list);
my $logfile= 'RAxML_logfile_' . $baselist;
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
my $topdir= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/RAxML_output/';
unless (-e $topdir){
	mkdir $topdir || die "Cannot create $topdir: $!\n";
}

for my $batch (sort {$a cmp $b} (keys %files)){
	for my $infile (sort {$a cmp $b} (keys %{$files{$batch}})){
		#parse the locus from full path to nucleotide alignment
		my $locus;
		my $basename= basename($infile);
		if ($basename=~ /^(HOG2_)(\d+?)(\.fa)/){
			$locus= $1 . $2;
		}
		else{
			print "Error parsing $locus from: $basename\n";
			die;
		}
		print "Processing: $locus\n";
	
		#if we've already aligned this locus...skip
		if (exists $done{$locus}){
			next;
		}
		else{
			#create a subdirectory for this batch...
			my $outdir= $topdir . $batch . '/';
			unless (-e $outdir){
				mkdir $outdir || die "Cannot create $outdir: $!\n";
			}
	
			#first, clean up any partially completed jobs before restarting (b/c RAxML will fail if output files already exist)
			my $out1= $outdir . "RAxML_bootstrap.$locus";
			if (-e $out1){
				unlink($out1) || die "Cannot delete $out1: $!\n";
			}
			my $out2= $outdir . "RAxML_info.$locus";
			if (-e $out2){
				unlink($out2) || die "Cannot delete $out2: $!\n";
			}
			my $out3= $outdir . "RAxML_bestTree.$locus";
			if (-e $out3){
				unlink($out3) || die "Cannot delete $out3: $!\n";
			}
			my $out4= $outdir . "RAxML_bipartitionsBranchLavels.$locus";
			if (-e $out4){
				unlink($out4) || die "Cannot delete $out4: $!\n";
			}
			my $out5= $outdir . "RAxML_bipartitions.$locus";
			if (-e $out5){
				unlink($out5) || die "Cannot delete $out5: $!\n";
			}

			#want to run partitioned by c1+c2/c3
			my $partitionfile= $outdir . $locus . '_partitions';
			#delete if already exists from a previous run...
			if (-e $partitionfile){
				unlink($partitionfile) || die "Cannot delete $partitionfile: $!\n";
			}
			#get the total alignment length...
			my $alignio_in= Bio::AlignIO-> new(-file=> $infile, -format=> 'fasta');
			my $aln_len;
			while (my $aln= $alignio_in-> next_aln()){
				$aln_len= $aln-> length();
			}
			#write partition file
			open(PARTITION, ">$partitionfile") || die "Cannot open $partitionfile for writing: $!\n";
			print PARTITION "DNA, gene1c1andc2 = 1-$aln_len\\3, 2-$aln_len\\3\n";
			print PARTITION "DNA, gene1c3 = 3-$aln_len\\3\n";
			close(PARTITION); 
		
			#specify substitution model
			#NB- RAxML default is GTR + G, but specify anyhow...
			my $model= 'GTRGAMMA';	
	
			#generate random numbers to use as the bootstrap seed and parsimony inference seed (for parsimony starting tree)
			#NB- you will have a record of these numbers in the .info output file from RAxML, but you could also output them to a 2nd logfile if desired
			my $range= 100000;
			my $min= 1;
			my $rand1= int(rand($range)) + $min;
			my $rand2= int(rand($range)) + $min;
		
			#command to run raxml
			#in this case, no outgroup, partitioned, 200 rapid bootstrap replicates and best ML tree, and with branch lengths calculated
			my $raxml= "raxmlHPC-SSE3 -f a -k -s $infile -p $rand1 -x $rand2 -# 200 -q $partitionfile -m $model -n $locus -w $outdir";

			#run as a system command
			system($raxml);

			#when finished, record in logfile so we'll skip redoing this locus if script stops/restarts
			print LOG "$locus done\n";
		}
	}
}
close(LOG);
