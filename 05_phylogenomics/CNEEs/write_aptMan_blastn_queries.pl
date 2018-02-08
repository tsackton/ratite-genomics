#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use File::Basename;

#create an output SeqIO stream to write retrieved sequences for each species
my $outfile1= 'aptHaa_blastn_queries.fasta';
my $seqio_out1= Bio::SeqIO-> new(-file=> ">$outfile1", -format=> 'fasta');

my $outfile2= 'aptOwe_blastn_queries.fasta';
my $seqio_out2= Bio::SeqIO-> new(-file=> ">$outfile2", -format=> 'fasta');

my $outfile3= 'aptRow_blastn_queries.fasta';
my $seqio_out3= Bio::SeqIO-> new(-file=> ">$outfile3", -format=> 'fasta');

my $output1= 0;
my $output2= 0;
my $output3= 0;
my $skipped= 0;
my $processed= 0;
for (my $i= 1; $i <= 15; $i++){
	print "Processing: batch$i\n";
	#glob alignment fastas in each batch
	my @glob= <CNEEs_unaligned/batch$i/*.fasta>;
	for my $filepath (sort {$a cmp $b} (@glob)){
		$processed++;
		my $basename= basename($filepath);
		my $locus;
		if ($basename=~ /^(.+)(\.fasta)$/){
			$locus= $1;
		}
		else{
			print "Error parsing locus from: $basename\n";
			die;
		}
	
		my $seqio_in= Bio::SeqIO-> new(-file=> $filepath, -format=> 'fasta');
		while (my $seqobj= $seqio_in-> next_seq()){
			my $id= $seqobj-> display_id();
			#skip if not a kiwi...
			unless ($id=~ /^(apt)/){
				next;
			}
			my $seq= $seqobj-> seq();
			#create a new seq object, with batch, then cnee as id
			my $newid= "batch$i\_$locus";
			my $seqobj2= Bio::Seq-> new(-id=> $newid, -seq=> $seq);
				
			#output to appropriate query seq file
			if ($id eq 'aptHaa'){
				$seqio_out1-> write_seq($seqobj2);
				$output1++;
			}
			elsif ($id eq 'aptOwe'){
				$seqio_out2-> write_seq($seqobj2);
				$output2++;
			}
			elsif ($id eq 'aptRow'){
				$seqio_out3-> write_seq($seqobj2);
				$output3++;
			}
			else{
				print "Error!  Unexpected taxon: $id\n";
				die;
			}	
		}
	}
}
print "Finished processing: $processed fastas\n";
print "Output: $output1 aptHaa seqs\n";
print "Output: $output2 aptOwe seqs\n";
print "Output: $output3 aptRow seqs\n\n";
