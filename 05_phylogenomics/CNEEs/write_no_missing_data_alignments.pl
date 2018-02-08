#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::SeqIO;
use Bio::Seq;

#top-level output directory
my $topdir= 'CNEEs_mafft_no-missing/';
unless (-e $topdir){
	mkdir($topdir) || die "Cannot create $topdir: $!\n";
}

for (my $batch= 1; $batch <= 15; $batch++){
	#glob all alignments
	my @glob= <CNEEs_mafft/batch$batch/*fasta>;
	my $scalarglob= scalar(@glob);
	print "Preparing to process: $scalarglob alignments in batch$batch\n";

	my $subdir= $topdir . 'batch' . $batch . '/';
	unless (-e $subdir){
		mkdir($subdir) || die "Cannot create $subdir: $!\n";
	}
	my $count= 0;
	my $output= 0;
	for my $filepath (sort {$a cmp $b} (@glob)){
		my $basename= basename($filepath);
		my $locus;
		if ($basename=~ /^(.+)(\.fasta)$/){
			$locus= $1;
		}
		else{
			print "Error parsing 'locus' from: $basename\n";
			die;
		}
		my $newaln1= Bio::SimpleAlign-> new();

		my $alnio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
		while(my $aln= $alnio_in-> next_aln()){
			if ($aln-> is_flush()){
				my $len= $aln-> length();
				my $numseqs= $aln-> num_sequences();
				#only proceed if there are no missing taxa...
				unless ($numseqs == 15){
					$count++;
					next;
				}
				else{
					#for each sequence...change all non-ACGT bases to - [gap]
					foreach my $seqobj ($aln-> each_seq()){
						my $id= $seqobj-> display_id();
						my $seq= $seqobj-> seq();
						$seq= uc($seq);
						$seq=~ s/[^ACGT]/-/g;
						my $locatableseq= Bio::LocatableSeq-> new(-id=> $id, -seq=> $seq);
						$newaln1-> add_seq($locatableseq);
					}
					#now, remove any columns containing 'gaps' (will remove all columns containing any gaps/missing/ambiguous sites)
					$newaln1= $newaln1-> remove_columns(['gaps']);
					my $newlen1= $newaln1-> length();
					#if we still have min. 200 bp...print out alignment
					if ($newlen1 >= 200){
						my $outfile= $subdir . $basename;
						my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
						foreach my $printobj ($newaln1-> each_seq()){
							$seqio_out-> write_seq($printobj);
						}
						$output++;
					}
					$count++;
					if ($count=~ /^(\d+)(0{2})$/){
						print "Processed: $count\n";
					}
				}
			}
			else{
				print "Error!  Alignment: $filepath is not flush!\n";
				die;
			}
		}
	}
	print "Finished processing: $count loci in batch$batch, Output: $output\n";
}
