#!/usr/bin/perl
use strict;
use warnings;

use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;

#list of all taxids
my @taxa= qw( anoDid aptHaa aptMan aptOwe aptRow casCas cryCin droNov eudEle galGal notPer rheAme rhePen strCam tinGut );
#initialize hash with each taxid as key, and empty string (to hold sequence) as value
my %seqs;
for my $spp (@taxa){
	$seqs{$spp}= '';
}

#go through each locus in a specified order...
my $count= 0;
my $totlen= 0;
#glob all input alignments...
my @glob= <introns_mafft/*.fasta>;
for my $filepath (sort {$a cmp $b} (@glob)){
	$count++;
	if ($count=~ /^(\d+)(0{2})/){
		print "Processed: $count\n";
	}
	my $alnio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
	while (my $aln= $alnio_in-> next_aln()){
		if ($aln-> is_flush()){
			my $aln_len= $aln-> length();
			for my $id (sort {$a cmp $b} (keys %seqs)){
				my $seqobj= $aln-> get_seq_by_id($id);
				if (defined($seqobj)){
					my $seq= $seqobj-> seq();
					$seqs{$id}= $seqs{$id} . $seq;
				}
				#if there's no sequence for this species, add all-gaps string...
				else{
					my $seq2= '-' x $aln_len;
					$seqs{$id}= $seqs{$id} . $seq2;
				}
			}
		}
		else{
			print "Error!  Alignment $filepath is not flush!\n";
				die;
		}
	}
}
#create an output alignment & then output in phylip format
my $newaln= Bio::SimpleAlign-> new();
for my $id2 (sort {$a cmp $b} (keys %seqs)){
	my $seq2= $seqs{$id2};
	$seq2= uc($seq2);
	my $locatable_seqobj= Bio::LocatableSeq-> new(-id=> $id2, -seq=> $seq2, -alphabet=> 'dna');
	$newaln-> add_seq($locatable_seqobj);
}
my $outfile2= 'introns_mafft_concat.phylip';
my $alnio_out= Bio::AlignIO-> new(-file=> ">$outfile2", -format=> 'phylip');
if ($newaln-> is_flush()){
	$totlen= $newaln-> length();
	$alnio_out-> write_aln($newaln);
}
else{
	print "Error!  Final alignment is not flush!\n";
	die;
}
print "Finished writing concatenated seqs. for $count loci\n";
print "Total concatenated alignment length: $totlen\n";
