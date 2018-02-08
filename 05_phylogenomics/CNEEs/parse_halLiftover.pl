#!/usr/bin/perl
use strict;
use warnings;

my $ref= 'galGal';

#array holding IDs of all taxa (except galGal reference)
my @taxa= qw( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );

for my $spp (@taxa){
	my %loci;
	my $index= 0;
	print "\nParsing: $spp liftover...\n";
	my $liftfile= "halLiftover/$ref-$spp\_cnees_halLiftover.psl";
	unless (-e $liftfile){
		print "Error!!! No file: $liftfile\n";
		die;
	}
	open(LIFT, "<$liftfile") || die "Cannot open $liftfile for reading: $!\n";
	while(<LIFT>){
		$index++;
		chomp($_);
		my @split= split(/\t/, $_);
		my ($locus, $scaff, $start, $end, $strand) = ($split[0], $split[14], $split[16], $split[17], $split[9]);
		
		$loci{$locus}{$index}{'scaffold'}= $scaff;
		$loci{$locus}{$index}{'start'}= $start;
		$loci{$locus}{$index}{'end'}= $end;
		if ($strand=~ /^([+-])([+-])$/){
			$loci{$locus}{$index}{'strand'}= $2;
			#expect all galGal strands to be +, but check...
			my $galstr= $1;
			unless ($galstr eq '+'){
				print "Error! $locus galGal strand is $galstr (expect +)\n";
				die;
			}
		}
		else{
			print "Error parsing strand from: $strand\n";
			die;
		}
	}
	close(LIFT);
	my %keep;
	#output all good liftover loci to .bed file for this species
	my $outfile= $spp . '_cnees_parsed_liftover.bed';
	open(OUTBED, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
	#only want to keep if there is no more than 1 seq/locus (e.g. keep only loci with a unique liftover region in target species)
	for my $key1 (sort {$a cmp $b} (keys %loci)){
		my $ckscalar= scalar(keys %{$loci{$key1}});
		if ($ckscalar == 1){
			$keep{$key1}= 'exists';
			for my $key2 (keys %{$loci{$key1}}){
				print OUTBED "$loci{$key1}{$key2}{'scaffold'}\t$loci{$key1}{$key2}{'start'}\t$loci{$key1}{$key2}{'end'}\tID=$key1\t0\t$loci{$key1}{$key2}{'strand'}\n";
			}
		}
	}
	close(OUTBED);
	my $scal1= scalar(keys %loci);
	my $scal2= scalar(keys %keep);
	print "Finished parsing: $scal1 loci for $spp...keep $scal2\n";	
}
