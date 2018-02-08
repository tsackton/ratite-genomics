#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;

#go through mafft alignment fastas subdirectories for each alignment batch
my $total= 0;
for (my $batchkey= 1; $batchkey <= 15; $batchkey++){
	my $checked= 0;
	
	print "\nWriting RAxML input files for: batch$batchkey\n";
	#create an output directory for this dataset/batch
	my $topdir= 'CNEEs_mafft/';
	unless (-e $topdir){
		mkdir($topdir) || die "Cannot create: $topdir: $!\n";
	}

	my $outdir= $topdir . 'batch' . $batchkey . '/';
	unless (-e $outdir){
		mkdir ($outdir) || die "Cannot create: $outdir: $!\n";
	}

	#glob all of the input fastas...
	my @glob= <CNEEs_mafft_temp/batch$batchkey/*.fasta>;
	my $scalarglob= scalar(@glob);

	for my $filepath (sort {$a cmp $b} (@glob)){
		my $basename= basename($filepath);
		#print "\nFile: $basename\n";
		my $locus;
		if ($basename=~ /^(.+)(_)(aln)(\.fasta)$/){
			$locus= $1;
		}
		else{
			print "Error parsing locus from: $basename\n";
			die;
		}

		#create an output SeqIO stream for writing
		my $outfile= $outdir . $locus . '.fasta';
		my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
		
		my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
		while (my $aln= $alignio_in-> next_aln()){
			if ($aln-> is_flush()){
				my %goodcol;
				#create a new alignment with all sequences in uppercase
				my $newaln= Bio::SimpleAlign-> new();
				#go through each seq. 
				foreach my $seqobj ($aln-> each_seq()){
					my $id= $seqobj-> display_id();
					my $seq= $seqobj-> seq();
					$seq= uc($seq);
					my $newseqobj= Bio::LocatableSeq-> new(-id=> $id, -seq=> $seq);
					$newaln-> add_seq($newseqobj);
					#split sequence & record any positions with bases other than - (gap) or N
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
				$checked++;
				$total++;
			}
			else{
				print "Error!  $locus alignment is not flush!\n";
				die;
			}
		}
	}
	print "\nProcessed: $checked fastas ($total total)\n\n";
}					
