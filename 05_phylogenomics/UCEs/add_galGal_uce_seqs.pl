#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use File::Copy;

#store all the galGal seqs. to add
my $addfile= 'galGal_uces.bed';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my %toadd;
while(<ADDSEQS>){
	chomp($_);
	my @split= split(/\t/, $_);
	my $idtemp= $split[3];
	my $query;
	if ($idtemp=~ /^(galGal,)(.+)$/){
		$query= $2;
	}
	else{
		print "Error parsing locus from: $idtemp\n";
		die;
	}
	my $hit= $split[0];
	my $hitstrand= $split[5];
	my $hitstart= $split[1];
	#add 1 to make 0-based bed start 1-based for fasta retrieval
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

#path to directory holding fastas
my $outdir= 'input_fastas_palaeos/';
my $processed= 0;
my $added= 0;

my @glob= <$outdir/*.fasta>;
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
	if (exists $toadd{$locname}){
		print "Adding galGal to: $locname\n";
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
		$added++;
	}
}
print "Finished processing: $processed fastas\n";
print "Added galGal to: $added\n\n";
