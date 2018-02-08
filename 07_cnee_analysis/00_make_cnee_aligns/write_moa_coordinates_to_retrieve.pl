#!/usr/bin/perl
use strict;
use warnings;

use Math::Round;

#append to existing logfile to record which seqs. are multiple/too long/too short/missing
open(LOG, ">moa_final_cnees_long_liftover_parsing_log") || die "Cannot open logfile for writing: $!\n";

open(LOG2, ">moa_cnee_parsing_log2") || die "Cannot open logfile2 for writing: $!\n";

#store lengths of all galGal reference cnees
my $dataset= 'final_cnees_long.bed';
my %cnees;
open(CNEES, "<$dataset") || die "Cannot open $dataset for reading: $1\n";
while(<CNEES>){
	chomp($_);
	my @spl= split(/\t/, $_);
	my ($pos1, $pos2, $name)= ($spl[1], $spl[2], $spl[3]);
	my $len= $pos2 - $pos1;
	#store by just numeric portion of name for easier sorting...
	my $numid;
	if ($name=~ /^([A-Za-z]+)(\d+)$/){
		$numid= $2;
	}
	else{
		print "Error parsing numerical portion of: $name\n";
		die;
	}
	$cnees{$numid}{'name'}= $name;
	$cnees{$numid}{'length'}= $len;
}
close(CNEES);
my $tot= scalar(keys %cnees);
print "\nFinished storing info. for: $tot cnees\n\n";

#read through the list of moa (anoDid) coverage (= % of galGal ref covered by non-N moa bases)
#store as too long if > 200% of galGal ref, and short if < 50%, only for cnees in 'final_cnees_long' dataset
#ALSO- if there's 0% coverage...store as 'no_liftover'
my $covfile= 'anoDid_cnee_coverage';
my %coverage;
open(COV, "<$covfile") || die "Cannot open $covfile for reading: $!\n";
while(<COV>){
	chomp($_);
	if ($_=~ /^(ID)/){
		#print "$_\n";
		next;
	}
	else{
		my @spl2= split(/\t/, $_);
		my ($cnee2, $cov2)= ($spl2[0], $spl2[6]);
		my $numid2;
		if ($cnee2=~ /^([A-Za-z]+)(\d+)$/){
			$numid2= $2;
		}
		if (exists $cnees{$numid2}){
			if ($cov2 > 200){
				$coverage{$cnee2}= 'target_region_too_long';
				#print "$_\n";
			}
			elsif (($cov2 > 0) && ($cov2 < 50)){
				$coverage{$cnee2}= 'target_region_short';
				#print "$_\n";
			}
			elsif ($cov2 == 0){
				$coverage{$cnee2}= 'no_liftover';
				#print "$_\n";
			}
			else{
				next;
			}
		}
	}
}
close(COV);
print "\nFinished storing info. about moa (anoDid) 'coverage' of galGal ref cnee length\n";
my $scalar2= scalar(keys %coverage);
print "$scalar2 cnees are too long/short in moa\n";

#store moa coordinates for all cnees in 'final_cnees_long' dataset
my $coordfile= 'moa-droNov_coordinates_to_retrieve';
my %coords;
open(COORDS, "<$coordfile") || die "Cannot open $coordfile for reading: $!\n";
while(<COORDS>){
	chomp($_);
	if ($_=~ /^(Locus)/){
		next;
	}
	else{
		my @spl3= split(/\t/, $_);
		my $cnee3= $spl3[0];
		my $numid3;
		if ($cnee3=~ /^([A-Za-z]+)(\d+)$/){
			$numid3= $2;
		}
		if (exists $cnees{$numid3}){
			$coords{$cnee3}= $_;
		}
	}
}
close(COORDS);
my $scalar3= scalar(keys %coords);
print "\nFinished storing coordinates to retrieve for: $scalar3 moa cnees\n";

my $checked= 0;
my $missing= 0;
my $toolong= 0;
my $tooshort= 0;
my $output= 0;

#output good liftover loci to file
my $outfile= 'moa-droNov_coordinates_to_retrieve_final';
open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
print OUT "Locus\tMoa_scaffold\tStart(1-based)\tEnd(1-based)\tStrand\n";

#go through ALL cnees (from galGal reference)
for my $keyb (sort {$a <=> $b} (keys %cnees)){
	$checked++;
	my $idb= $cnees{$keyb}{'name'};	

	#check if we have this cnee in target species
	unless (exists $coords{$idb}){
		print LOG "$idb\tanoDid\tno_liftover\n";
		$missing++;
		next;
	}
	else{
		my $printline= $coords{$idb};
		#check the start/end coordinates...
		my @spl4= split(/\t/, $printline);
		my $start4= $spl4[2];
		my $end4= $spl4[3];
		if ($start4 >= $end4){
			print LOG "$idb\tanoDid\tno_liftover\n";
			$missing++;
			next;
		}
		#if we flagged this cnee as too long/short...
		if (exists $coverage{$idb}){
			my $comment= $coverage{$idb};
			#if too long...omit
			if ($comment eq 'target_region_too_long'){
				print LOG "$idb\tanoDid\t$comment\n";
				$toolong++;
			}
			#if too short, print comment to log, but don't omit (so print out coords)
			elsif ($comment eq 'target_region_short'){
				print LOG "$idb\tanoDid\t$comment\n";
				$tooshort++;
				print OUT "$printline\n";
				$output++;
			}
			#if coverage was 0...output as missing (no liftover)
			elsif ($comment eq 'no_liftover'){
				print LOG "$idb\tanoDid\t$comment\n";
				$missing++;
			}
			else{
				print "Error: unexpected comment for $idb ($comment)\n";
				die;
			}
		}
		#if not flagged for coverage, just print out
		else{
			print OUT "$printline\n";
			$output++;
		}
	}
}
close(OUT);
print "Finished processing: $checked loci for anoDid\n";
print "No liftover region: $missing (omitted)\n";
print "Overly long target region (> 2X galGal): $toolong (omitted)\n";
print "Short target region (< 50% of chicken): $tooshort (retained)\n";
print "Output coordinates for: $output anoDid cnees\n\n";

print LOG2 "Finished processing: $checked loci for anoDid\n";
print LOG2 "No liftover region: $missing (omitted)\n";
print LOG2 "Overly long target region (> 2X galGal): $toolong (omitted)\n";
print LOG2 "Short target region (< 50% of chicken): $tooshort (retained)\n";
print LOG2 "Output coordinates for: $output anoDid cnees\n\n";

close(LOG);
close(LOG2);
