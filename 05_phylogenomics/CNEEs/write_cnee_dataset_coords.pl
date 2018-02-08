#!/usr/bin/perl
use strict;
use warnings;

my @taxa= qw( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );
my $scalarspp= scalar(@taxa);

#store the 'good' liftovers for each taxon...
my %loci;
for my $spp (@taxa){
	print "\nProcessing: $spp\n";
	my $bedfile= "halLiftover/$spp\_cnees_parsed_liftover.bed";
	unless (-e $bedfile){
		print "Error!  No file: $bedfile\n";
		die;
	}
	open(BEDIN, "<$bedfile") || die "Cannot open $bedfile for reading: $!\n";
	my $count= 0;
	while(<BEDIN>){
		chomp($_);
		my @split= split(/\t/, $_);
		my $len= $split[2] - $split[1];
		#set a minimum seq. length cutoff of 250 bp
		unless ($len >= 250){
			next;
		}
		if ($split[3]=~ /^(ID=)(.+)$/){
			$loci{$2}{$spp}= 'exists';
			$count++;
		}
		else{
			print "Error parsing id: $split[3]\n";
			die;
		}
	}
	close(BEDIN);
	print "Stored info for: $count loci\n";
}

#keep only loci that have sequence for ALL palaeognaths
my $goodloc= 0;
my %tokeep;
for my $key1 (sort {$a cmp $b} (keys %loci)){
	my $scalartax= scalar(keys %{$loci{$key1}});
	if ($scalartax == $scalarspp){
		$goodloc++;
		$tokeep{$key1}= 'exists';
	}
}
print "Good loci (0 missing spp., min. length 250 bp): $goodloc\n";
my $scalar2= scalar(keys %tokeep);
print "Check: $scalar2\n";

#print a list of all of loci to use
open(OUTLIST, ">candidate_cnees_list") || die "Cannot open outfile for writing: $!\n";
for my $cle (sort {$a cmp $b} (keys %tokeep)){
	print OUTLIST "$cle\n";
}
close(OUTLIST);

#output .bed lines for loci we want to keep...
for my $spp2 (@taxa){
	my $cnt= 0;
	my $bedfile2= "halLiftover/$spp2\_cnees_parsed_liftover.bed";
	my $outfile= "halLiftover/$spp2\_cnees_to_use.bed";
	open(BEDIN2, "<$bedfile2") || die "Cannot open $bedfile2 for reading: $!\n";
	open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
	
	while(<BEDIN2>){
		chomp($_);
		my @split2= split(/\t/, $_);
		if ($split2[3]=~ /^(ID=)(.+)$/){
			my $id2= $2;
			if (exists $tokeep{$id2}){
				print OUT "$_\n";
				$cnt++;
			}
		}
		else{
			print "Error parsing: $split2[3]\n";
			die;
		}
	}
	close(BEDIN2);
	close(OUT);
	print "Output: $cnt loci for $spp2\n";
}
