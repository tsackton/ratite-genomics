#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#create a fastadb of the moa-droNov genome seq.
my $fasta= '/n/regal/edwards_lab/acloutier/moa_mapping/moa-droNov_bowtie_remapping/moa-droNov_bowtie_remapping_consensus_genome_final.fasta';
print "Creating fastadb of: $fasta\n";
my $fastadb= Bio::DB::Fasta-> new($fasta);

#create an output fasta to write moa cnees
my $seqio_out= Bio::SeqIO-> new(-file=> ">anoDid_cnees.fa", -format=> 'fasta');

#go through the coordinates for all the moa (anoDid) seqs. to add
my $addfile= 'moa-droNov_coordinates_to_retrieve';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my $processed= 0;
my $added= 0;
my $skipped= 0;
while(<ADDSEQS>){
	chomp($_);
	if ($_=~ /^(Locus)/){
		next;
	}
	else{
		$processed++;
		my @split= split(/\t/, $_);
		my ($query, $hit, $hitstart, $hitend, $hitstrand)= ($split[0], $split[1], $split[2], $split[3], $split[4]);
		
		#if start is greater than or equal to end...skip
		if ($hitstart >= $hitend){
			print "Skipping: $query (start: $hitstart, end: $hitend)\n";
			$skipped++;
			next;
		}

		#get moa scaffold
		my $seqobj= $fastadb-> get_Seq_by_id($hit);
		unless (defined($seqobj)){
			print "Skipping: $query (no seqobj for $hit)\n";
			$skipped++;
			next;
		}
		#get the subseq (NB- 'hitstart' is already 1-based coordinate)
		my $subseq= $seqobj-> subseq($hitstart,$hitend);
		
		#create a new seqobj
		my $newid= 'anoDid_' . $query;
		my $new_seqobj= Bio::Seq-> new(-id=> $newid, -seq=> $subseq, -alphabet=> 'dna');
		#reverse complement if needed...
		my $final_seqobj;
		if ($hitstrand eq '+'){
			$final_seqobj= $new_seqobj;
		}
		elsif ($hitstrand eq '-'){
			$final_seqobj= $new_seqobj-> revcom();
		}
		else{
			print "Error!  Unexpected strand: $hitstrand\n";
		}

		$seqio_out-> write_seq($final_seqobj);
		$added++;
	}
}
close(ADDSEQS);
print "Finished processing: $processed CNEEs\n";
print "Output: $added moa CNEE seqs, skipped $skipped\n\n";
