#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;

#get current batch from command line
my $batch= shift(@ARGV);

#open output file to hold partition info (RAxML-style partition file)
open(OUT, ">$batch\_partitions") || die "Cannot open outfile to write partitions: $!\n";

#store cnee_batches
my %batches;
my $listfile= 'concat_lists/list' . $batch;
open(BATCHES, "<$listfile") || die"Cannot open list of batches for reading:$!\n";
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
		if ($count=~ /^(\d+)(0{3})$/){
			print "Processed: $count\n";
		}
		#path to alignment for this locus
		my $alnfile= "input_fastas_allspecies_cnees_aligned_no-galGal-gaps/$batchkey/$lockey\.fasta";
		my $alnio_in= Bio::AlignIO-> new(-file=> $alnfile, -format=> 'fasta');
		while (my $aln= $alnio_in-> next_aln()){
			if ($aln-> is_flush()){
				my $aln_len= $aln-> length();
				my $start= $totlen + 1;
				#my $bedstart= $totlen;
				$totlen+= $aln_len;
				print OUT "DNA, $lockey = $start-$totlen\n";
			}
			else{
				print "Error!  Alignment $lockey is not flush!\n";
				die;
			}
		}
	}
}
close(OUT);
print "Finished writing partitions for $count loci\n";
print "Total concatenated alignment length: $totlen\n";
