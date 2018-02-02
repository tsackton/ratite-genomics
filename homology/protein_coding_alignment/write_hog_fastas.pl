#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;

#file obtained from Tim Sackton: /n/edwards_lab/Users/tsackton/ratites/homology/post_hmm/new_hog_list.txt
my $infile= 'new_hog_list.txt';
open(INFILE, "<$infile") || die "Cannot open $infile for reading: $!\n";

my %hogs;

while(<INFILE>){
	chomp($_);
	my @split= split(/\t/, $_);
	my ($hog, $prot, $spp)= ($split[0], $split[1], $split[3]);
	$hogs{$hog}{$spp}{$prot}= 'exists';
}
close(INFILE);
my $scalar= scalar(keys %hogs);
print "\nStored: $scalar hogs\n";
my %toget;
my $sppcount= 0;
my $count= 0;
#only want to process hogs that have at least 4 spp.
for my $key1 (sort {$a cmp $b} (keys %hogs)){
	$count++;
	my $numspp= scalar(keys %{$hogs{$key1}});
	if ($numspp < 4){
		$sppcount++;
	}
	else{
		$toget{$key1}= 'exists';
	}
}
print "\nFinished processing: $count hogs:\n";
print "$sppcount have < 4 spp. total (skip)\n";
my $scalarkeep= scalar(keys %toget);
print "Preparing to output fastas for: $scalarkeep hogs\n\n";

my $fastafile= 'all_proteins.fa';
print "Creating fastadb...\n";
my $fastadb= Bio::DB::Fasta-> new($fastafile);

#create a toplevel directory to hold output fastas
my $topdir= 'fastas/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}

my $total= 0;
my $output= 0;
my $index= 0;
my $subdir;
for my $keya (sort {$a cmp $b} (keys %hogs)){
	if (exists $toget{$keya}){
		#create a new output subdirectory if needed
		if ($output == 0){
			$index++;
			$subdir= $topdir . 'batch' . $index . '/';
			unless (-e $subdir){
				mkdir $subdir || die "Cannot create $subdir: $!\n";
			}
		}
		#create an output SeqIO stream
		my $outfile= $subdir . $keya . '.fa';
		my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
		for my $keyb (sort {$a cmp $b} (keys %{$hogs{$keya}})){
			for my $keyc (sort {$a cmp $b} (keys %{$hogs{$keya}{$keyb}})){
				my $seqobj= $fastadb-> get_Seq_by_id($keyc);
				#change the id to have spp name, then protid
				my $id= $keyb . '_' . $keyc;
				my $seq= $seqobj-> seq();
				my $new_seqobj= Bio::Seq-> new(-id=> $id, -seq=> $seq);
				$seqio_out-> write_seq($new_seqobj);
			}
		}
		$output++;
		#write progress statement so we aren't wondering what's happening...
		if ($output =~ /^(\d+)(0{2})$/){
			print "Output: $output\n";
		}
		$total++;
		if ($output == 4500){
			$output= 0;
		}
	}
}
print "Finished outputting: $total hogs to $index subdirectories\n\n";
