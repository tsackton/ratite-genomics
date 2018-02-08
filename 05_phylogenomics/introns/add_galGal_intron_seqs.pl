#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use File::Copy;

#store all the galGal 'nonoverlapping' intron coordinates to add
my $addfile= 'galGal4_nonoverlapping_introns.bed';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my %toadd;
while(<ADDSEQS>){
	chomp($_);
	my @split= split(/\t/, $_);
	
	my $idtemp= $split[3];
	my $query;
	if ($idtemp=~ /^(CDS=)(.+?)(,)(.+)(intron=)(\d+)$/){
		$query= $2 . '_intron' . $6;
	}
	else{
		print "Error parsing locus from: $idtemp\n";
		die;
	}
	my $hit= $split[0];
	my $hitstrand= $split[5];
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

#directory holding fastas (we'll append galGal to these)
my $outdir= 'introns_unaligned_all/';

my $processed= 0;
for (my $i= 1; $i<= 34; $i++){
	print "Processing: batch$i\n";
	#specify subdirectory for this batch...
	my $subdir= $outdir . 'batch' . $i . '/';
	unless (-e $subdir){
		mkdir $subdir || die "Cannot create $subdir: $!\n";
	}

	#glob all the input fastas
	my @glob= <introns_unaligned_all/batch$i/*.fasta>;
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

		#if this is a locus that has an galGal seq to add (all should)...
		if (exists $toadd{$locname}){
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
			my $seqio_out= Bio::SeqIO-> new(-file=> ">>$filepath", -format=> 'fasta');
			$seqio_out-> write_seq($final_seqobj);
		}
		else{
			print "Error: no stored info for $locname\n";
		}
	}
}
print "Finished processing: $processed fastas\n";
