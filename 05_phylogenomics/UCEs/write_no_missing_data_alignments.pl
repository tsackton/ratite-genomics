#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::SeqIO;
use Bio::Seq;

#store any completed jobs
my %donejobs;
my $logfile= 'uces_no-missing_data_alignments_logfile';
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		my @spl= split(/\t/, $_);
		$donejobs{$spl[0]}= 'exists';
	}
	close(DONE);
	open(LOG, ">>$logfile") || die "Cannot open $logfile for writing: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

#glob all input mafft alignments
my @glob= <UCEs_mafft/*fasta>;
my $scalarglob= scalar(@glob);
print "Preparing to process: $scalarglob alignments\n";

#top-level output directory
my $topdir= 'UCEs_mafft_no-missing/';
unless (-e $topdir){
	mkdir($topdir) || die "Cannot create $topdir: $!\n";
}

my $count= 0;
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
	if (exists $donejobs{$locus}){
		$count++;
		next;
	}
	my $newaln1= Bio::SimpleAlign-> new();

	my $alnio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
	while(my $aln= $alnio_in-> next_aln()){
		if ($aln-> is_flush()){
			my $len= $aln-> length();
			my $numseqs= $aln-> num_sequences();
			#skip if there are any missing taxa
			unless ($numseqs == 15){
				print LOG "$locus\n";
				next;
			}
			else{
				foreach my $seqobj ($aln-> each_seq()){
					my $id= $seqobj-> display_id();
					my $seq= $seqobj-> seq();
					$seq= uc($seq);
					#overwrite any non-ACGT bases with gap characters
					$seq=~ s/[^ACGT]/-/g;
					my $locatableseq= Bio::LocatableSeq-> new(-id=> $id, -seq=> $seq, -alphabet=> 'dna');
					$newaln1-> add_seq($locatableseq);
				}
				#remove any columns containing 'gaps' (will remove all columns containing any gaps/missing/ambiguous sites)
				$newaln1= $newaln1-> remove_columns(['gaps']);
				my $newlen1= $newaln1-> length();
				#min threshold length 200 bp after removal of columns with missing data
				if ($newlen1 >= 200){
					my $outfile= $topdir . $basename;
					my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
					foreach my $printobj ($newaln1-> each_seq()){
						$seqio_out-> write_seq($printobj);
					}
				}
				print LOG "$locus\n";			
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
print "Finished processing: $count loci\n";
