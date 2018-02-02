#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

my $infile= 'good_PAML_HOG_protein_transcript_info';
open(INFILE, "<$infile") || die "Cannot open $infile for reading: $!\n";

my %hogs;
my %hogstemp;
my %hogs2;
my %hogs3;
my $hogcount= 0;
my $batch= 0;
my $totseqs= 0;

while(<INFILE>){
	chomp($_);
	if ($_=~ /^(HOG\t)/){
		next;
	}
	$totseqs++;
	my @split= split(/\t/, $_);
	my ($hog, $spp, $transcript)= ($split[0], $split[1], $split[4]);
	#append spp_ to start of transcript ID so will match entry in cds fasta
	$transcript= $spp . '_' . $transcript;
	$hogs{$spp}{$hog}{$transcript}= 'exists';
	$hogstemp{$hog}= 'exists';	
}
close(INFILE);
#assign each HOG to a batch...
for my $hkey (sort {$a cmp $b} (keys %hogstemp)){
	if ($hogcount == 0){
		$batch++;
		$hogs2{$hkey}= $batch;
		$hogs3{$batch}{$hkey}= 'exists';
		$hogcount++;
	}
	elsif ($hogcount <= 999){
		$hogs2{$hkey}= $batch;
		$hogs3{$batch}{$hkey}= 'exists';
		$hogcount++;
		if ($hogcount == 1000){
			$hogcount= 0;
			next;
		}
	}
	else{
		print "Error: unexpected 'hogcount' value: $hogcount\n";
		die;
	}
}

my $scalar= scalar(keys %hogs2);
my $scalar2= scalar(keys %hogs);
print "\nStored info for: $totseqs seqs. for $scalar2 species in $scalar hogs in $batch batches\n";
#output a list of which batch each HOG is in
open(BATCHLIST, ">HOG_cds_batches") || die "Cannot open file to write list of HOG batches: $!\n";
print BATCHLIST "Batch\tHOG\n";
for my $k1 (sort {$a <=> $b} (keys %hogs3)){
	for my $k2 (sort {$a cmp $b} (keys %{$hogs3{$k1}})){
		print BATCHLIST "batch$k1\t$k2\n";
	}
}
close(BATCHLIST);

#create a toplevel directory to hold output fastas
my $topdir= 'PAML_input_PRANK_fastas/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}

my $output= 0;
my $totoutput= 0;
my $subdir;
for my $keya (sort {$a cmp $b} (keys %hogs)){
	print "Processing: $keya\n";
	
	#create a fastadb of the cds seqs
	my $fastafile= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/cds_fastas/' . $keya . '_cds.fa'; 
	my $fastadb= Bio::DB::Fasta-> new($fastafile);

	for my $keyb (sort {$a cmp $b} (keys %{$hogs{$keya}})){
		my $batchdir;
		if (exists $hogs2{$keyb}){
			$batchdir= $hogs2{$keyb};
		}	
		else{
			print "Error!  No batch info. stored for: $keyb\n";
			die;
		}
			
		#create a new output subdirectory if needed
		$subdir= $topdir . 'batch' . $batchdir . '/';
		unless (-e $subdir){
			mkdir $subdir || die "Cannot create $subdir: $!\n";
		}

		#create an output SeqIO stream
		my $outfile= $subdir . $keyb . '.fa';
		#open for writing/appending
		my $seqio_out;
		if (-e $outfile){
			$seqio_out= Bio::SeqIO-> new(-file=> ">>$outfile", -format=> 'fasta');
		}
		else{
			$seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
		}
		#output each of the sequences for this species in this HOG group
		for my $keyc (sort {$a cmp $b} (keys %{$hogs{$keya}{$keyb}})){
			my $seqobj= $fastadb-> get_Seq_by_id($keyc);
			unless(defined($seqobj)){
				print "Error: $keyb no $keya seq $keyc...SKIPPED\n";
			}
			else{
				$seqio_out-> write_seq($seqobj);
				$output++;
				$totoutput++;
			}
		}
	}
	print "Output: $output seqs for $keya\n";
	$output= 0;
}
print "Output: $totoutput seqs total\n\n";
