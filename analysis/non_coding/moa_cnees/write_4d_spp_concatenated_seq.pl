#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;

#get the species from the command line
my $spp= shift(@ARGV);

#fastadb of genome
my $fasta;
if ($spp eq 'anoDid'){
	$fasta= '/n/regal/edwards_lab/acloutier/moa_mapping/moa-droNov_bowtie_remapping/moa-droNov_bowtie_remapping_consensus_genome_final.fasta';
}
else{
	$fasta= '/n/regal/edwards_lab/acloutier/refGenomes/' . $spp . '.fa';
}
my $fastadb= Bio::DB::Fasta-> new($fasta);

#string to hold concatenated 4d site info
my $seq= '';

#path to .bed file holding 4d sites
my $bedfile;
my $output= 0;
my $retrieved= 0;
my $missing= 0;
my $multiple= 0;
if ($spp eq 'galGal'){
	$bedfile= 'galGal4_4d_bed6_withid.bed';
}
else{
	$bedfile= "$spp\.4dparsed.bed";
}
open(BED, "<$bedfile") || die "Cannot open $bedfile for reading: $!\n";
while(<BED>){
	chomp($_);
	#if 'missing', append a gap to string
	if ($_=~ /^(missing)/){
		$seq= $seq . '-';
		$missing++;
		$output++;
	}
	#if 'multiple', append an N
	elsif ($_=~ /^(multiple)/){
		$seq= $seq . 'N';
		$multiple++;
		$output++;
	}
	#otherwise...retrieve scaffold
	else{
		my @split= split(/\t/, $_);
		my ($scaff, $start, $end, $strand)= ($split[0], $split[1], $split[2], $split[5]);
		my $seqobj= $fastadb-> get_Seq_by_id($scaff);
		unless(defined($seqobj)){
			print "Couldn't retrieve $scaff for $spp (output as 'missing')\n";
			$seq= $seq . '-';
			$missing++;
			$output++;
			next;
		}
		else{
			#increment 0-based bed start
			$start++;
			my $subseq= $seqobj-> subseq($start,$end);
			my $outbase;
			if ($strand eq '+'){
				$outbase= $subseq;
			}
			elsif ($strand eq '-'){
				my $new_seqobj= Bio::Seq-> new(-id=> $spp, -seq=> $subseq, -alphabet=> 'dna');
				$new_seqobj= $new_seqobj-> revcom();
				$outbase= $new_seqobj-> seq();
			}
			else{
				print "Error!  Unexpected strand: $strand\n";
				die;
			}
			$seq= $seq . $outbase;
			$retrieved++;
			$output++;
		}
	}
}
close(BED);
$seq= uc($seq);
my $final_seqobj= Bio::Seq-> new(-id=> $spp, -seq=> $seq, -alphabet=> 'dna');
my $fasta_out= $spp . '_4d_concatenated.fasta';
my $seqio_out= Bio::SeqIO-> new(-file=> ">$fasta_out", -format=> 'fasta');
$seqio_out-> write_seq($final_seqobj);
my $seqlen= length($seq);
print "Output $spp length: $seqlen\n";
print "Retrieved: $retrieved bases\n";
print "Missing: $missing (output as -)\n";
print "Multiple: $multiple (output as N)\n";
my $totoutput= $retrieved + $missing + $multiple;
unless (($totoutput == $seqlen) && ($totoutput == $output)){
	print "Error! Output $output bases, but seq length is $seqlen & (retrieved+missing+multiple) is $totoutput\n";
}
