#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#store all cnee_batches
my %batches;
open(BATCHES, 'cnee_alignment_batches') || die"Cannot open list of batches for reading:$!\n";
while(<BATCHES>){
	chomp($_);
	if ($_=~ /^(Batch)/){
		next;
	}
	else{
		my @spl= split(/\t/, $_);
		#store keyed by batch, then locus
		$batches{$spl[0]}{$spl[1]}= 'exists';
	}
}
close(BATCHES);

#list of all taxa (including chicken reference [galGal], but not moa [anoDid])
my @taxa= qw( allMis allSin anaPla anoCar aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla galGal gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut );

#get the batch number from command line (for running split into multiple jobs)
my $batch= shift(@ARGV);

my $topdir= 'input_fastas_allspecies_cnees/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}
my $outdir= $topdir . $batch . '/';
unless (-e $outdir){
	mkdir($outdir) || die "Cannot create $outdir: $!\n";
}

#store any already completed jobs (from stop/restart on serial requeue partition)
my $logfile= '/n/regal/edwards_lab/acloutier/moa_cnees/write_fastas_log_' . $batch;
my %donejobs;
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(.+?)(\s)(done)$/){
			$donejobs{$1}= 'exists';
		}
	}
	close(DONE);
	#reopen file for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}
#delete any output fastas that aren't stored in %donejobs
print "\nCleaning up files from Odyssey stop/restart...\n";
my @glob= </n/regal/edwards_lab/acloutier/moa_cnees/$outdir/*.fasta>;
for my $filepath (sort {$a cmp $b} (@glob)){
	my $basename= basename($filepath);
	my $rootname;
	if ($basename=~ /^(.+)(\.)(fasta)$/){
		$rootname= $1;
	}
	else{
		print "Error parsing 'root' name from: $basename\n";
		die;
	}
	unless (exists $donejobs{$rootname}){
		unlink($filepath);
	}
}

#for each species...
for my $spp (@taxa){
	print "Processing: $spp\n";
	#fastadb of genome
	my $fasta= '/n/regal/edwards_lab/acloutier/refGenomes/' . $spp . '.fa';
	my $fastadb= Bio::DB::Fasta-> new($fasta);

	my $coordfile;
	if ($spp eq 'galGal'){
		$coordfile= 'final_cnees_long.bed';
	}
	else{
		$coordfile= $spp . '_cnees_parsed_liftover.bed';
	}
	#for each cnee present in this species...
	my %coords;
	open(COORDS, "<$coordfile") || die "Cannot open $coordfile for reading: $!\n";
	while(<COORDS>){
		chomp($_);
		my @split2= split(/\t/, $_);
		my ($scaff, $start, $end, $cnee, $strand)= ($split2[0], $split2[1], $split2[2], $split2[3], $split2[5]);
		#increment start to be 1-based
		$start++;
		#set galGal strand to + for all loci
		if ($spp eq 'galGal'){
			$strand= '+';
		}
		$coords{$cnee}{'scaffold'}= $scaff;
		$coords{$cnee}{'start'}= $start;
		$coords{$cnee}{'end'}= $end;
		$coords{$cnee}{'strand'}= $strand;
	}
	close(COORDS);
	my $numloci= scalar(keys %coords);
	print "Stored: $numloci\n";

	#for each cnee in this batch...
	for my $cneekey (sort {$a cmp $b} (keys %{$batches{$batch}})){
		#skip if we've already processed locus...
		if (exists $donejobs{$cneekey}){
			next;
		}
		my $outfile= $outdir . $cneekey . '.fasta';
			
		#if file already exists, open for appending; otherwise open for writing
		my $seqio_out;
		if (-e $outfile){
			$seqio_out= Bio::SeqIO-> new(-file=> ">>$outfile", -format=> 'fasta');
			#print "Appending to: $outfile\n";
		}
		else{
			$seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
			#print "Writing to: $outfile\n";
		}
		if (exists $coords{$cneekey}){
			my $scaff2= $coords{$cneekey}{'scaffold'};
			my $start2= $coords{$cneekey}{'start'};
			my $end2= $coords{$cneekey}{'end'};
			my $strand2= $coords{$cneekey}{'strand'};

			my $seqobj= $fastadb-> get_Seq_by_id($scaff2);
			unless(defined($seqobj)){
				print "Error!  No scaffold $scaff2 for $spp ($batch, $cneekey)\n";
				next;
			}
			my $subseq= $seqobj-> subseq($start2,$end2);
			my $new_seqobj= Bio::Seq-> new(-id=> $spp, -seq=> $subseq, -alphabet=> 'dna');
			#reverse complement if on - strand
			if ($strand2 eq '-'){
				$new_seqobj= $new_seqobj-> revcom();
			}
			$seqio_out-> write_seq($new_seqobj);
		}
		#if this is the last species (tinGut), record that this locus is done
		if ($spp eq 'tinGut'){
			print LOG "$cneekey done\n";
		}
	}
	close(COORDS);
}
