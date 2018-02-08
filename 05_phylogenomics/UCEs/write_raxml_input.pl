#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;

my $outdir= 'UCEs_mafft_raxml_input/';
unless (-e $outdir){
	mkdir($outdir) || die "Cannot create: $outdir: $!\n";
}

#glob all of the input fastas...
my @glob= <UCEs_mafft/*.fasta>;
my $scalarglob= scalar(@glob);

for my $filepath (sort {$a cmp $b} (@glob)){
	if (-s $filepath){
		my $basename= basename($filepath);
		my $locus;
		if ($basename=~ /^(.+)(_)(aln)(\.fasta)$/){
			$locus= $1;
		}
		else{
			print "Error parsing locus from: $basename\n";
			die;
		}
		my $outfile= $outdir . $locus . '.fasta';
		my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');

		my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
		while (my $aln= $alignio_in-> next_aln()){
			if ($aln-> is_flush()){
				my %goodcol;

				my $newaln= Bio::SimpleAlign-> new();
				foreach my $seqobj ($aln-> each_seq()){
					my $id= $seqobj-> display_id();
					my $seq= $seqobj-> seq();
					$seq= uc($seq);
					my $newseqobj= Bio::LocatableSeq-> new(-id=> $id, -seq=> $seq);
					$newaln-> add_seq($newseqobj);
					#now, split sequence & record any positions with bases other than - (gap) or N
					my @split= split(//, $seq);
					my $index= 0;
					for my $base (@split){
						if (($base ne '-') && ($base ne 'N')){
							#store 'good' position as 0-based value
							$goodcol{$index}= 'exists';
						}
						$index++;
					}
				}
				my $alnlen= $newaln-> length();
				my $scalargood= scalar(keys %goodcol);

				#go through each alignment column beginning at END
				for (my $i= $alnlen; $i>= 1; $i--){
					#get the 0-based equivalent...
					my $col= $i - 1;
					#if we had any good bases...keep
					if (exists $goodcol{$col}){
						next;
					}
					#otherwise, remove...
					else{
						$newaln= $newaln-> remove_columns([$col,$col]);
					}
				}
				#double-check that final alignment length is as expected
				my $newlen= $newaln-> length();
				unless ($newlen == $scalargood){
					print "Error! $locus: orig. length= $alnlen, 'good' columns (to keep)= $scalargood, output length= $newlen\n";
					die;
				}
				#output alignment as a SeqIO stream...
				foreach my $seqobj2 ($newaln-> each_seq()){
					$seqio_out-> write_seq($seqobj2);
					my $id2= $seqobj2-> display_id();
				}
			}
			else{
				print "Error!  $locus alignment is not flush!\n";
				die;
			}
		}
	}
}				
