#!/usr/bin/perl
use strict;
use warnings;

use POSIX;

my @taxa= qw ( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );
my $numtaxa= scalar(@taxa);

my %loci;
my %strcam_len;
for my $spp (@taxa){
	print "\nProcessing: $spp...\n";
	my $bedfile= 'halLiftover/' . $spp . '_nonoverlapping_introns.bed';
	open(BED, "<$bedfile") || die "Cannot open $bedfile for reading: $!\n";
	while(<BED>){
		chomp($_);
		my @split= split(/\t/, $_);
		my $id= $split[3];
		my $cds;
		my $intron;
		
		if ($id=~ /^(CDS=)(.+?)(,intron=)(\d+)$/){
			$cds= $2;
			$intron= $4;
		}
		else{
			print "Error parsing: $id\n";
			die;
		}
		#only store as a good intron if tot. length is at least 100 bp
		my $chlen= $split[2] - $split[1];
		if ($chlen >= 100){
			$loci{$cds}{$intron}{$spp}= 'exists';
		}
		#if current spp. is strCam, also store the intron length
		if ($spp eq 'strCam'){
			my $len= $split[2] - $split[1];
			$strcam_len{$cds}{$intron}= $len;
		}
	}
	close(BED);
}
my $scalar1= scalar(keys %loci);
print "\nFinshed storing intron coordinates for: $scalar1 CDS\n";

#retain introns with seq for all species
my $count= 0;
my %tokeep;
my %cdslen;
for my $cdskey (sort {$a cmp $b} (keys %loci)){
	for my $intronkey (sort {$a <=> $b} (keys %{$loci{$cdskey}})){
		my $scalar_taxa= scalar(keys %{$loci{$cdskey}{$intronkey}});
		if ($scalar_taxa == $numtaxa){
			#get the length of this 'good' intron in strCam...
			my $checklen= $strcam_len{$cdskey}{$intronkey};
			if (exists $cdslen{$cdskey}){
				$cdslen{$cdskey}+= $checklen;
			}
			else{
				$cdslen{$cdskey}= $checklen;
			}
			$count++;
			$tokeep{$cdskey}{$intronkey}= 'exists';
		}
	}
}
my $scalar2= scalar(keys %tokeep);
print "\nRetained: $count introns from $scalar2 CDS that have intron(s) with seqs. for all taxa\n\n";

#now, filter out cases where total CDS length in strCam is < 200 bp
print "Filtering CDS with total strCam length < 200 bp\n\n";
my %tokeep2;
my $count2= 0;
for my $cdskey2 (sort {$a cmp $b} (keys %tokeep)){
	my $totlen= $cdslen{$cdskey2};
	if ($totlen >= 200){
		for my $intronkey2 (keys %{$tokeep{$cdskey2}}){
			$tokeep2{$cdskey2}{$intronkey2}= 'exists';
			$count2++;
		}
	}
}
my $scalar2b= scalar(keys %tokeep2);
print "\nRetained: $count2 introns from $scalar2b CDS that meet min. length requirement\n\n";

#last, in the case of alternative transcripts, pick the longest one (based on strCam), or else pick one randomly
print "Choosing 1 CDS among alternative transcripts for loci...\n\n";

my $galfile= 'intput_annotations/galGal4_nonoverlapping_introns.bed';
open(GAL, "<$galfile") || die "Cannot open $galfile for reading: $!\n";
my %transcripts;
while(<GAL>){
	chomp($_);
	my @split2= split(/\t/, $_);
	my $id2= $split2[3];
	my $cds2;
	my $gene2;
	if ($id2=~ /^(CDS=)(.+?)(,gene=)(.+?)(,)/){
		$cds2= $2;
		$gene2= $4;
		$transcripts{$gene2}{$cds2}= 'exists';
	}
	else{
		print "Error parsing galGal annotation: $id2\n";
	}
}
close(GAL);
my $scal= scalar(keys %transcripts);
print "\nFinished storing alternative transcript info. for $scal galGal genes\n";

my %tokeep3;
my $count3= 0;
for my $genekey3 (sort {$a cmp $b} (keys %transcripts)){
	my $numcds= scalar(keys %{$transcripts{$genekey3}});
	#if there's a single CDS for this gene, and it's one we have stored 'to keep'...
	if ($numcds == 1){
		for my $cdskey3 (keys %{$transcripts{$genekey3}}){
			if (exists $tokeep2{$cdskey3}){
				for my $intronkey3 (keys %{$tokeep2{$cdskey3}}){
					$tokeep3{$cdskey3}{$intronkey3}= 'exists';
					$count3++;
				}
			}
		}
	}
	#if there are multiple transcripts for this gene...
	else{
		#store the total 'good' length for each CDS we actually have stored in a temp hash
		my %temphash;
		for my $cdskey4 (keys %{$transcripts{$genekey3}}){
			if (exists $cdslen{$cdskey4}){
				my $ln= $cdslen{$cdskey4};
				$temphash{$ln}{$cdskey4}= 'exists';
			}
		}
		my $numkept= scalar(keys %temphash);
		
		#if we didn't have any 'good' CDS anyhow...
		if ($numkept == 0){
			next;
		}
		#otherwise, for the longest 'good' CDS length...
		else{
			my @lens= sort {$a <=> $b} (keys %temphash);
			my $longest= pop(@lens);
	
			my $numkept2= scalar(keys %{$temphash{$longest}});
			#if there's only a single 'good' CDS...use it
			if ($numkept2 == 1){
				for my $usekey2 (keys %{$temphash{$longest}}){
					#loop through the stored introns from this CDS...
					for my $usekey3 (keys %{$tokeep2{$usekey2}}){
						$tokeep3{$usekey2}{$usekey3}= 'exists';
						$count3++;
					}
				}
			}
			#if there are multiple 'good' CDS with this same longest strCam length, pick one at random
			else{
				my @cdsids= keys(%{$temphash{$longest}});
				my $range= scalar(@cdsids)- 1;
				my $random_number= int(rand($range));
				my $useid= $cdsids[$random_number];
				#store the introns to keep
				for my $keepkey1 (keys %{$tokeep2{$useid}}){
					$tokeep3{$useid}{$keepkey1}= 'exists';
					$count3++;
				}
			}
		}
	}
}
my $scalar3= scalar(keys %tokeep3);
print "\nRetained: $count3 introns from $scalar3 CDS after keeping only 1 alternative transcript/gene\n\n";

open(OUT, ">candidate_intron_list") || die "Cannot open outfile for writing: $!\n";
print OUT "CDS\tintron\n";
for my $printkey1 (sort {$a cmp $b} (keys %tokeep3)){
	for my $printkey2 (sort {$a <=> $b} (keys %{$tokeep3{$printkey1}})){
		print OUT "$printkey1\t$printkey2\n";
	}
}
close(OUT);
