#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use File::Copy;

#store all the moa (anoDid) seqs. to add
my $addfile= '/n/regal/edwards_lab/acloutier/moa_cnees/cnee_coordinates/moa-droNov_coordinates_to_retrieve_final';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my %toadd;
while(<ADDSEQS>){
	chomp($_);
	if ($_=~ /^(Locus)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		my $query= $split[0];
		my $hit= $split[1];
		my $hitstrand= $split[4];
		my $hitstart= $split[2];
		my $hitend= $split[3];

		$toadd{$query}{'hit'}= $hit;
		$toadd{$query}{'strand'}= $hitstrand;
		$toadd{$query}{'start'}= $hitstart;
		$toadd{$query}{'end'}= $hitend;
	}
}
close(ADDSEQS);
my $scalaradd= scalar(keys %toadd);
print "\nFinished storing: $scalaradd seqs. to add...\n\n";

#create a fastadb of the moa-droNov genome seq.
my $fasta= '/n/regal/edwards_lab/acloutier/moa_mapping/moa-droNov_bowtie_remapping/moa-droNov_bowtie_remapping_consensus_genome_final.fasta';
print "Creating fastadb of: $fasta\n";
my $fastadb= Bio::DB::Fasta-> new($fasta);

#path to directory holding fastas
my $outdir= 'input_fastas_allspecies_cnees/';

#glob all of the input fastas in each batch...
my $processed= 0;
my $added= 0;
for (my $i= 1; $i<= 57; $i++){
	print "Processing: batch$i\n";
	#specify subdirectory for this batch...
	my $subdir= $outdir . 'batch' . $i . '/';
	unless (-e $subdir){
		mkdir $subdir || die "Cannot create $subdir: $!\n";
	}

	#glob all the input fastas
	my @glob= </n/regal/edwards_lab/acloutier/moa_cnees/input_fastas_allspecies_cnees/batch$i/*.fasta>;
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
		#if this is a locus that has an anoDid seq to add...
		if (exists $toadd{$locname}){
		
			my $getseq= $toadd{$locname}{'hit'};
			my $getstrand= $toadd{$locname}{'strand'};
			my $getstart= $toadd{$locname}{'start'};
			my $getend= $toadd{$locname}{'end'};

			my $seqobj= $fastadb-> get_Seq_by_id($getseq);
			unless(defined($seqobj)){
				print "No scaffold: $getseq\n";
				next;
			}
			#get the subseq
			my $subseq= $seqobj-> subseq($getstart,$getend);
			
			#create a new seqobj
			my $new_seqobj= Bio::Seq-> new(-id=> 'anoDid', -seq=> $subseq, -alphabet=> 'dna');
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

			#create a SeqIO output stream to APPEND anoDid seq. to fasta
			my $seqio_out= Bio::SeqIO-> new(-file=> ">>$filepath", -format=> 'fasta');
			$seqio_out-> write_seq($final_seqobj);
			$added++;
		}
	}
}
print "Finished processing: $processed fastas\n";
print "Added moa (anoDid) to: $added\n\n";
