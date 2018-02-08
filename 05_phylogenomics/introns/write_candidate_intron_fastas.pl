#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

my @taxa= qw ( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );

#get list of all introns for a given 'batch'
my $infile= shift(@ARGV);
my $basename= basename($infile);

#parse the 'batch' number from input intron list
my $batch;
if ($basename=~ /^(intron_list)(\d+)$/){
	$batch= $2;
}
else{
	print "Error parsing 'batch' number from: $infile\n";
	die;
}
my $topdir= 'introns_unaligned_all/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}
my $outdir= $topdir . 'batch' . $batch . '/';
unless (-e $outdir){
	mkdir($outdir) || die "Cannot create $outdir: $!\n";
}

#store any already completed jobs
my $logfile= 'write_fastas_log' . $batch;
my %donejobs;
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(.+?)(\s)(.+?)(\s)(done)$/){
			$donejobs{$1}{$3}= 'exists';
		}
	}
	close(DONE);
	#reopen file for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

#go through the introns we're going to process, but skip any we've already completed...
my %introns;
open(INFILE, "<$infile") || die "Cannot open $infile for reading: $!\n";
while(<INFILE>){
	chomp($_);
	if ($_=~ /^(cds)/){
		my @split= split(/\t/, $_);
		if (exists $donejobs{$split[0]}{$split[1]}){
			next;
		}
		$introns{$split[0]}{$split[1]}= 'exists';
	}
}
close(INFILE);

my %coords;
for my $spp (@taxa){
	print "Processing: $spp\n";
	#create fastadb for each species
	my $fasta= 'genomes/' . $spp . '.fa';
	my $fastadb= Bio::DB::Fasta-> new($fasta);

	my $coordfile= 'halLiftover/' . $spp . '_nonoverlapping_introns.bed';
	open(COORDS, "<$coordfile") || die "Cannot open $coordfile for reading: $!\n";
	while(<COORDS>){
		chomp($_);
		my @split2= split(/\t/, $_);
		my $id2= $split2[3];
		my $cds2;
		my $intron2;
		if ($id2=~ /^(CDS=)(.+?)(,intron=)(\d+)$/){
			$cds2= $2;
			$intron2= $4;
		}
		else{
			print "Error parsing: $id2\n";
			die;
		}
		#if this is a CDS/intron we're going to use...
		if (exists $introns{$cds2}{$intron2}){
			my $outfile= $outdir . $cds2 . '_intron' . $intron2 . '.fasta';
			#if file already exists, open for appending; otherwise open for writing
			my $seqio_out;
			if (-e $outfile){
				$seqio_out= Bio::SeqIO-> new(-file=> ">>$outfile", -format=> 'fasta');
			}
			else{
				$seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
			}

			#retrieve the intron sequence from fasta
			my $scaff= $split2[0];
			my $start= $split2[1];
			#increment start to be 1-based
			$start++;
			my $end= $split2[2];
			my $strand= $split2[5];

			my $seqobj= $fastadb-> get_Seq_by_id($scaff);
			my $subseq= $seqobj-> subseq($start,$end);
			my $new_seqobj= Bio::Seq-> new(-id=> $spp, -seq=> $subseq);
			#reverse complement if on - strand
			if ($strand eq '-'){
				$new_seqobj= $new_seqobj-> revcom();
			}
			$seqio_out-> write_seq($new_seqobj);
			#if this is the last species (tinGut), record that this locus is done
			if ($spp eq 'tinGut'){
				print LOG "$cds2 $intron2 done\n";
			}
		}
		
	}
	close(COORDS);
}
