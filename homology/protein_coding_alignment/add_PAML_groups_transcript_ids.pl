#!/usr/bin/perl
use strict;
use warnings;

#store all of the protid & transcript ids from Tim's all_info file, keyed first by spp id
#Tim's file accessed from: /n/edwards_lab/Users/tsackton/ratites/homology/all_info.summary
my $infofile= 'all_info.summary';
open(INFO, "<$infofile") || die "Cannot open $infofile for reading: $!\n";
my %idinfo;
my %allspp;
while(<INFO>){
	chomp($_);
	if ($_=~ /^(#species)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		my ($spp, $prot, $trans)= ($split[0], $split[9], $split[6]);
		$allspp{$spp}= 'exists';
		#store if protein id is not <NA>, and has the expected NCBI accession.version format OR expected maker format (for team ratite genomes)
		#e.g. also ignore protein IDs like 'Ig lambda chain V-IV region Bau-like'
		if (($prot ne '<NA>') && (($prot=~ /^([A-Z]{2})(_)(\d+)(\.)(\d+)$/) || ($prot=~ /^([A-Z]{4})(\d{8})(-)([A-Z]{2})$/))){
			#NCBI GFFs have protein accession.version, but Tim's HOG list just had accession...
			my $shortprot;
			if ($prot=~ /^(.+?)(\.)(\d+)$/){
				$shortprot= $1;
			}
			elsif ($prot=~ /^([A-Z]{4})(\d{8})(-)([A-Z]{2})$/){
				$shortprot= $prot;
			}
			else{
				print "Error parsing protein ACCESSION (minus version) from: $prot\n";
				die;
			}
			$idinfo{$spp}{$shortprot}{$trans}= $prot;
		}
	}
}
close(INFO);
#check all is as expected...
my $scalarspp= scalar(keys %idinfo);
print "Stored protein & transcript IDs for: $scalarspp species\n";
for my $k1 (sort {$a cmp $b} (keys %allspp)){
	unless (exists $idinfo{$k1}){
		print "No stored info. for: $k1\n";
	}
}
for my $key1 (sort {$a cmp $b} (keys %idinfo)){
	for my $key2 (sort {$a cmp $b} (keys %{$idinfo{$key1}})){
		my $scalarids= scalar(keys %{$idinfo{$key1}{$key2}});
		unless ($scalarids == 1){
			print "\nSpp: $key1, protein: $key2 has $scalarids transcript ids:\n";
			for my $key3 (sort {$a cmp $b} (keys %{$idinfo{$key1}{$key2}})){
				print "$key3\n";
				#die;
			}
		}
	}
}

#now, read in the good PAML HOG groups & output both versions of protein ID, and the corresponding transcript ID
open(OUT, ">good_PAML_HOG_protein_transcript_info") || die "Cannot open outfile for writing: $!\n";
print OUT "HOG\tSpecies\tHOG_protein_id\tNCBI_protein_id\ttranscript_id\n";

my $infile= '/home/alison/Desktop/oma_duplicate_genes/good_PAML_HOG_protein_info';
open(INFILE, "<$infile") || die "Cannot open $infile for reading: $!\n";
my $lines= 0;
while(<INFILE>){
	chomp($_);
	if ($_=~ /^(HOG\t)/){
		next;
	}
	else{
		$lines++;
		my @split2= split(/\t/, $_);
		my ($hog2, $spp2, $prot2)= ($split2[0], $split2[1], $split2[2]);
		if (exists $idinfo{$spp2}{$prot2}){
			for my $trans2 (keys %{$idinfo{$spp2}{$prot2}}){
				my $longprot= $idinfo{$spp2}{$prot2}{$trans2};
				print OUT "$hog2\t$spp2\t$prot2\t$longprot\t$trans2\n";
			}
		}
		else{
			print "Error!  No info. for: $spp2 $prot2\n";
			die;
		}
	}
}
close(INFILE);
close(OUT);
print "Added transcript info to: $lines protein ids\n\n";
