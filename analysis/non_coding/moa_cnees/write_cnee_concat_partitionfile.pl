#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;

#open output file to hold partition info (RAxML-style partition file)
open(OUT, ">allspecies_cnee_concat_partitions") || die "Cannot open outfile to write partitions: $!\n";

#2nd output file to hold partition info in .bed format
open(BEDOUT, ">allspecies_cnee_concat_partitions.bed") || die "Cannot open .bed outfile to write partitions: $!\n";

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

#go through each locus in each batch...in same order as for all species concatenation
my $count= 0;
my $totlen= 0;
for my $batchkey (sort {$a cmp $b} (keys %batches)){
	for my $lockey (sort {$a cmp $b} (keys %{$batches{$batchkey}})){
		$count++;
		#path to alignment for this locus
		my $alnfile= "input_fastas_allspecies_cnees_aligned_no-galGal-gaps/$batchkey/$lockey\.fasta";
		my $alnio_in= Bio::AlignIO-> new(-file=> $alnfile, -format=> 'fasta');
		while (my $aln= $alnio_in-> next_aln()){
			if ($aln-> is_flush()){
				my $aln_len= $aln-> length();
				my $start= $totlen + 1;
				my $bedstart= $totlen;
				$totlen+= $aln_len;
				print OUT "DNA, $lockey = $start-$totlen\n";
				print BEDOUT "$lockey\t$bedstart\t$totlen\t$lockey\n";
			}
			else{
				print "Error!  Alignment $lockey is not flush!\n";
				die;
			}
		}
	}
}
close(OUT);
close(BEDOUT);
print "Finished writing partitions for $count loci\n";
print "Total concatenated alignment length: $totlen\n";
