#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

my @taxa= qw ( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );

#CNEEs split into 'batches' to avoid directories with too many files...
my $infile= shift(@ARGV);
my $basename= basename($infile);

#parse the 'batch' number from input cnee list
my $batch;
if ($basename=~ /^(cnee_list)(\d+)$/){
	$batch= $2;
}
else{
	print "Error parsing 'batch' number from: $infile\n";
	die;
}
my $topdir= 'CNEEs_unaligned/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}
my $outdir= $topdir . 'batch' . $batch . '/';
unless (-e $outdir){
	mkdir($outdir) || die "Cannot create $outdir: $!\n";
}

my %cnees;
open(INFILE, "<$infile") || die "Cannot open $infile for reading: $!\n";
while(<INFILE>){
	chomp($_);
	$cnees{$_}= 'exists';
}
close(INFILE);

#store all good cnee coords for all species
my %coords;
for my $spp (@taxa){
	print "Processing: $spp\n";
	my $fasta= 'genomes/' . $spp . '.fa';
	my $fastadb= Bio::DB::Fasta-> new($fasta);

	my $coordfile= 'halLiftover/' . $spp . '_cnees_to_use.bed';
	open(COORDS, "<$coordfile") || die "Cannot open $coordfile for reading: $!\n";
	while(<COORDS>){
		chomp($_);
		my @split2= split(/\t/, $_);
		my $id2= $split2[3];
		my $cnee;
		if ($id2=~ /^(ID=)(.+)$/){
			$cnee= $2;
		}
		
		if (exists $cnees{$cnee}){
			my $outfile= $outdir . $cnee . '.fasta';
			#if file already exists, open for appending; otherwise open for writing
			my $seqio_out;
			if (-e $outfile){
				$seqio_out= Bio::SeqIO-> new(-file=> ">>$outfile", -format=> 'fasta');
			}
			else{
				$seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
			}

			#retrieve the cnee sequence from fasta
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
		}
		
	}
	close(COORDS);
}
