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
my $outfile= 'introns_mafft_concat_partitions';
open(OUT, ">$outfile") || die "Cannot open outfile to write partitions: $!\n";

#go through each locus in same order as for concatenated alignment
my $count= 0;
my $totlen= 0;
my @glob= <introns_mafft/*.fasta>;
for my $filepath (sort {$a cmp $b} (@glob)){
	$count++;
	if ($count=~ /^(\d+)(0{2})/){
		print "Processed: $count\n";
	}
	my $basename= basename($filepath);
	my $locus;
	if ($basename=~ /^(.+)(\.fasta)$/){
		$locus= $1;
	}
	else{
		print "Error parsing locus from: $basename\n";
		die;
	}
	my $alnio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
	while (my $aln= $alnio_in-> next_aln()){
		if ($aln-> is_flush()){
			my $aln_len= $aln-> length();
			my $start= $totlen + 1;
			$totlen+= $aln_len;
			print OUT "DNA, $locus = $start-$totlen\n";
		}
		else{
			print "Error!  Alignment $locus is not flush!\n";
			die;
		}
	}
}
close(OUT);
print "Finished writing partitions for $count loci\n";
print "Total concatenated alignment length: $totlen\n";
