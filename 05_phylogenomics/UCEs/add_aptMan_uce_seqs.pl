#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;
use File::Copy;

#store all the aptMan seqs. to add
my $addfile= 'aptMan_seqs_to_retrieve';
open (ADDSEQS, "<" . $addfile) || die "Cannot open $addfile for reading: $!\n";
my %toadd;
while(<ADDSEQS>){
	chomp($_);
	if ($_=~ /^(Query)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		my $query= $split[0];
		my $hit= $split[2];
		my $hitstrand= $split[3];
		my $hitstart= $split[4];
		my $hitend= $split[5];

		$toadd{$query}{'hit'}= $hit;
		$toadd{$query}{'strand'}= $hitstrand;
		$toadd{$query}{'start'}= $hitstart;
		$toadd{$query}{'end'}= $hitend;
	}
}
close(ADDSEQS);
my $scalaradd= scalar(keys %toadd);
print "\nFinished storing: $scalaradd seqs. to add...\n\n";

#create a fastadb of the aptMan genome seq.
my $fasta= 'genomes/aptMan.fa';
print "Creating fastadb of: $fasta\n";
my $fastadb= Bio::DB::Fasta-> new($fasta);

#directory holding fastas
my $outdir= 'UCEs_unaligned/';

my $processed= 0;
my $added= 0;

#glob all the input fastas
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
	#if this is a locus that has an aptMan seq to add...
	if (exists $toadd{$locname}){
		print "Adding aptMan to: $locname\n";
		my $getseq= $toadd{$locname}{'hit'};
		my $getstrand= $toadd{$locname}{'strand'};
		my $getstart= $toadd{$locname}{'start'};
		my $getend= $toadd{$locname}{'end'};

		my $seqobj= $fastadb-> get_Seq_by_id($getseq);
		#get the subseq
		my $subseq= $seqobj-> subseq($getstart,$getend);
		#create a new seqobj
		my $new_seqobj= Bio::Seq-> new(-id=> 'aptMan', -seq=> $subseq);
		#reverse complement if needed...
		my $final_seqobj;
		if ($getstrand == 1){
			$final_seqobj= $new_seqobj;
		}
		elsif ($getstrand == -1){
			$final_seqobj= $new_seqobj-> revcom();
		}
		else{
			print "Error!  Unexpected strand: $getstrand\n";
		}

		#create a SeqIO output stream to APPEND aptMan seq. to fasta
		my $seqio_out= Bio::SeqIO-> new(-file=> ">>$filepath", -format=> 'fasta');
		$seqio_out-> write_seq($final_seqobj);
		$added++;
	}
}
print "Finished processing: $processed fastas\n";
print "Added aptMan to: $added\n\n";
