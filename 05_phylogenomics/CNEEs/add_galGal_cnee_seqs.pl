#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use File::Copy;

#store all the galGal candidate cnee coordinates
my $addfile= 'galGal_candidate_cnees.bed';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my %toadd;
while(<ADDSEQS>){
	chomp($_);
	my @split= split(/\t/, $_);
	
	my $query= $split[3];
	my $hit= $split[0];
	my $hitstrand= '+';
	my $hitstart= $split[1];
	#add 1 to make 0-based bed start 1-based for fasta retrieval
	$hitstart++;
	my $hitend= $split[2];

	$toadd{$query}{'hit'}= $hit;
	$toadd{$query}{'strand'}= $hitstrand;
	$toadd{$query}{'start'}= $hitstart;
	$toadd{$query}{'end'}= $hitend;

}
close(ADDSEQS);
my $scalaradd= scalar(keys %toadd);
print "\nFinished storing: $scalaradd seqs. to add...\n\n";

#create a fastadb of the galGal genome seq.
my $fasta= 'genomes/galGal.fa';
print "Creating fastadb of: $fasta\n";
my $fastadb= Bio::DB::Fasta-> new($fasta);

#path to directory holding input palaeo. fastas
my $indir= 'input_fastas/';
#create an output directory
my $outdir= 'CNEEs_unaligned/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}

#glob all of the input fastas in each batch...
my $totadded= 0;
my $processed= 0;
my $couldadd= 0;
my $added= 0;
for (my $i= 1; $i<= 15; $i++){
	print "Processing: batch$i\n";
	#specify subdirectory for this batch...
	my $subdir= $outdir . 'batch' . $i . '/';
	unless (-e $subdir){
		mkdir $subdir || die "Cannot create $subdir: $!\n";
	}

	#glob all the input fastas
	my @glob= <input_fastas/batch$i/*.fasta>;
	my $scalarglob= scalar(@glob);
	print "Preparing to process: $scalarglob fastas\n";

	for my $filepath (sort {$a cmp $b} (@glob)){
		$processed++;
		my $basename= basename($filepath);
		my $locname;
		if ($basename=~ /^(.+?)(\.)(fasta)$/){
			$locname= $1;
		}
		else{
			print "Error parsing locus name from: $basename\n";
			die;
		}
		#copy the fasta of palaeo seqs.
		my $outfile= $subdir . $basename;
		copy($filepath,$outfile) || die "Copy failed for: $locname: $!\n";

		#if this is a locus that has an galGal seq to add (all should)...
		if (exists $toadd{$locname}){
			$couldadd++;
		
			my $getseq= $toadd{$locname}{'hit'};
			my $getstrand= $toadd{$locname}{'strand'};
			my $getstart= $toadd{$locname}{'start'};
			my $getend= $toadd{$locname}{'end'};

			my $seqobj= $fastadb-> get_Seq_by_id($getseq);
			#get the subseq
			my $subseq= $seqobj-> subseq($getstart,$getend);
			
			#create a new seqobj
			my $new_seqobj= Bio::Seq-> new(-id=> 'galGal', -seq=> $subseq);
			#reverse complement if needed...
			my $final_seqobj;
			if ($getstrand eq '+'){
				$final_seqobj= $new_seqobj;
			}
			elsif ($getstrand eq '-'){
				$final_seqobj= $new_seqobj-> revcom();
			}
			else{
				print "Error!  Unexpected strand: $getstrand\n";
			}

			#create a SeqIO output stream to APPEND galGal seq. to fasta
			my $seqio_out= Bio::SeqIO-> new(-file=> ">>$outfile", -format=> 'fasta');
			$seqio_out-> write_seq($final_seqobj);
			$added++;
			$totadded++;
		}
		else{
			print "Error: no stored info for $locname\n";
			die;
		}
	}
}
print "Finished processing: $processed fastas\n";
print "Considered adding galGal to: $couldadd loci\n";
print "Added galGal to: $added\n\n";
