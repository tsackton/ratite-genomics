#!/usr/bin/perl
use strict;
use warnings;

use Math::Round;

#create a logfile to record which seqs. are multiple/too long/too short/missing
open(LOG, ">final_cnees_long_liftover_parsing_log") || die "Cannot open logfile for writing: $!\n";
print LOG "Locus\tTaxon\tNote\n";

open(LOG2, ">cnee_parsing_log2") || die "Cannot open logfile2 for writing: $!\n";

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

#list of all taxa (except reference chicken galGal & moa anoDid)
my @taxa= qw( allMis allSin anaPla anoCar aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut );

for my $spp (@taxa){
	print "\nParsing liftover for: $spp\n";
	my %loci;
	my $index= 0;
	my $checked= 0;
	my $missing= 0;
	my $multiple= 0;
	my $toolong= 0;
	my $tooshort= 0;
	my $output= 0;

	my $liftfile= "$spp\.psl";		
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
			#expect all reference galGal strands to be +, but check...
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
	my $liftovers= scalar(keys %loci);

	#output good liftover loci to .bed file
	my $outfile= "$spp\_cnees_parsed_liftover.bed";
	open(OUTBED, ">$outfile") || die "Cannot open $outfile for writing: $!\n";

	#go through ALL cnees (from galGal reference)
	for my $keyb (sort {$a <=> $b} (keys %cnees)){
		$checked++;
		my $idb= $cnees{$keyb}{'name'};
		my $lenb= $cnees{$keyb}{'length'};	

		#check if we have this cnee in target species
		#if no liftover...
		unless (exists $loci{$idb}){
			print LOG "$idb\t$spp\tno_liftover\n";
			$missing++;
			next;
		}
		else{
			#check how many liftover targets we had for this cnee
			my $cknum= scalar(keys %{$loci{$idb}});
			#if we have a single liftover region...
			if ($cknum == 1){
				#check length relative to galGal
				for my $keyc (keys %{$loci{$idb}}){
					my $lenc= $loci{$idb}{$keyc}{'end'} - $loci{$idb}{$keyc}{'start'};
					my $cklen= round(($lenc/$lenb)*100);
					#if target region is more than 200% galGal length...omit
					if ($cklen > 200){
						print LOG "$idb\t$spp\ttarget_region_too_long\n";
						$toolong++;
						next;
					}
					#if target region is less than 50% of galGal, print note to log, but don't omit
					elsif ($cklen < 50){
						print LOG "$idb\t$spp\ttarget_region_short\n";
						$tooshort++;
						print OUTBED "$loci{$idb}{$keyc}{'scaffold'}\t$loci{$idb}{$keyc}{'start'}\t$loci{$idb}{$keyc}{'end'}\t$idb\t0\t$loci{$idb}{$keyc}{'strand'}\n";
						$output++;
					}
					#if target region is a reasonable length...just output bed coordinates
					else{
						print OUTBED "$loci{$idb}{$keyc}{'scaffold'}\t$loci{$idb}{$keyc}{'start'}\t$loci{$idb}{$keyc}{'end'}\t$idb\t0\t$loci{$idb}{$keyc}{'strand'}\n";
						$output++;
					}
				}
			}
			#if we had multiple liftover regions...
			else{
				print LOG "$idb\t$spp\tmultiple_liftover_regions\n";
				$multiple++;
				next;
			}
		}
	}
	close(OUTBED);

	print LOG2 "Finished processing: $checked loci for $spp\n";
	print LOG2 "No liftover region: $missing (omitted)\n";
	print LOG2 "Have cnee liftover for: $liftovers loci\n";
	print LOG2 "Multiple liftover regions: $multiple (omitted)\n";
	print LOG2 "Overly long target region (> 2X galGal): $toolong (omitted)\n";
	print LOG2 "Short target region (< 50% of chicken): $tooshort (retained)\n";
	print LOG2 "Output bed coordinates for: $output cnees for $spp\n\n";
}
close(LOG);
close(LOG2);
