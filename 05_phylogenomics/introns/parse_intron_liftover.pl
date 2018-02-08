#!/usr/bin/perl
use strict;
use warnings;

#get species from command line (e.g. droNov)
my $spp= shift(@ARGV);

#store the annotated boundaries for galGal4 'nonoverlapping' introns
my $galfile= 'input_annotations/galGal4_nonoverlapping_introns.bed';
open(GALREF, "<$galfile") || die "Cannot open $galfile for reading: $!\n";
my %gal;

my $count= 0;
while(<GALREF>){
	chomp($_);
	my @split= split(/\t/, $_);
	my $id= $split[3];
	my $cds;
	my $intron;
	if ($id=~ /^(CDS=)(.+?)(,)(.+?)(intron=)(\d+)$/){
		$cds= $2;
		$intron= $6;
		$count++;
	}
	else{
		print "Error parsing: $id\n";
		die;
	}

	#store info for scaffold, strand, start & end
	$gal{$cds}{$intron}{'scaffold'}= $split[0];
	$gal{$cds}{$intron}{'start'}= $split[1];
	$gal{$cds}{$intron}{'end'}= $split[2];
	$gal{$cds}{$intron}{'strand'}= $split[5];
}
close(GALREF);
my $scalar1= scalar(keys %gal);
print "\nFinished storing: $count introns from $scalar1 reference galGal CDS\n\n";

#now, store the expected target intron coordinates, based on matches to adjacent exons
my $targetfile= 'halLiftover/' . $spp . '_expected_intron_coords';
my %target;

open(TARGET, "<$targetfile") || die "Cannot open $targetfile for reading: $!\n";
my $count2= 0;

while(<TARGET>){
	chomp($_);
	if ($_=~ /^(cds)/){
		my @split2= split(/\t/, $_);
		my $tstart= $split2[4];
		my $tend= $split2[5];
		my $sflush= 'yes';
		my $eflush= 'yes';

		#if there's a > or < at start of intron start/end coord, exon match was not flush...
		if ($tstart=~ /^([^\d])(\d+)$/){
			$sflush= 'no';
			$tstart= $2;
		}
		if ($tend=~ /^([^\d])(\d+)$/){
			$eflush= 'no';
			$tend= $2;
		}

		$target{$split2[0]}{$split2[1]}{'scaffold'}= $split2[2];
		$target{$split2[0]}{$split2[1]}{'start'}= $tstart;
		$target{$split2[0]}{$split2[1]}{'end'}= $tend;
		$target{$split2[0]}{$split2[1]}{'strand'}= $split2[3];
		$target{$split2[0]}{$split2[1]}{'start_flush'}= $sflush;
		$target{$split2[0]}{$split2[1]}{'end_flush'}= $eflush;
		$count2++;
	}
}
close(TARGET);
my $scalar2= scalar(keys %target);
print "Finished storing: $count2 introns from $scalar2 target CDS\n\n";

#now, go through the halLiftover of introns between galGal & target spp
my $liftover= 'halLiftover/' . $spp . '_nonoverlapping_introns_liftover.psl';
open(LIFT, "<$liftover") || die "Cannot open $liftover for reading: $!\n";
#for sanity's sake, loop through a first time to check how many liftover 'matches' occur in the expected region
#and only use cases with a single match to update coords
my %regions;
my $count3= 0;
while(<LIFT>){
	chomp($_);
	my @split4= split(/\t/, $_);
	my $id4= $split4[0];
	my $cds4;
	my $intron4;
	if ($id4=~ /^(CDS=)(.+?)(,)(.+?)(intron=)(\d+)$/){
		$cds4= $2;
		$intron4= $6;
	}
	else{
		print "Error parsing: $id4\n";
		die;
	}
	
	#if we have this intron stored as one to use for target
	if (exists $target{$cds4}{$intron4}){
		#if the liftover region falls within that expected...
		my $i_scaff4= $split4[14];
		my $i_start4= $split4[16];
		my $i_end4= $split4[17];
		#parse the target strand
		my @strands4= split('', $split4[9]);
		my $tstrand4= pop(@strands4);

		my $t_start4= $target{$cds4}{$intron4}{'start'};
		my $t_end4= $target{$cds4}{$intron4}{'end'};

		#if current liftover region is on expected scaffold & strand in target AND within the expected region...
		if (($i_scaff4 eq $target{$cds4}{$intron4}{'scaffold'}) && ($tstrand4 eq $target{$cds4}{$intron4}{'strand'}) && ($i_start4 >= $t_start4) && ($i_end4 <= $t_end4)){
					
			$count3++;
			#intialize/increment count of how many 'matches' there are within this region...
			if (exists $regions{$cds4}{$intron4}){
				$regions{$cds4}{$intron4}++;
			}
			else{
				$regions{$cds4}{$intron4}= 1;
			}
		}
	}
}
close(LIFT);
my $scalar3= scalar(keys %regions);
print "\nFinished storing # liftover matches for: $count3 introns from $scalar3 CDS\n";

#loop through file a 2nd time to refine expected intron coordinates
my $checked= 0;
my $updated_start= 0;
my $updated_end= 0;
open(LIFT2, "<$liftover") || die "Cannot open $liftover for reading: $!\n";
while(<LIFT2>){
	chomp($_);
	my @split3= split(/\t/, $_);
	my $id3= $split3[0];
	my $cds3;
	my $intron3;
	if ($id3=~ /^(CDS=)(.+?)(,)(.+?)(intron=)(\d+)$/){
		$cds3= $2;
		$intron3= $6;
	}
	else{
		print "Error parsing: $id3\n";
		die;
	}
	
	#if we have this intron stored as one to use for target
	if (exists $target{$cds3}{$intron3}){
		#if the liftover region falls within that expected...
		my $i_scaff= $split3[14];
		my $i_start= $split3[16];
		my $i_end= $split3[17];
		#parse the query/target strand
		my @strands= split('', $split3[9]);
		my $qstrand= shift(@strands);
		my $tstrand= shift(@strands);

		#check if the query (galGal) match region is flush to the annotated intron start/end
		my $qstart= $split3[12];
		my $qend= $split3[13];
		my $q_sflush;
		my $q_eflush;
		if ($qstart == $gal{$cds3}{$intron3}{'start'}){
			$q_sflush= 'yes';
		}
		else{
			$q_sflush= 'no';
		}
		if ($qend == $gal{$cds3}{$intron3}{'end'}){
			$q_eflush= 'yes';
		}
		else{
			$q_eflush= 'no';
		}
		my $t_start= $target{$cds3}{$intron3}{'start'};
		my $t_end= $target{$cds3}{$intron3}{'end'};
		
		#if current liftover region is on expected scaffold & strand in target AND within the expected region...
		if (($i_scaff eq $target{$cds3}{$intron3}{'scaffold'}) && ($tstrand eq $target{$cds3}{$intron3}{'strand'}) && ($i_start >= $t_start) && ($i_end <= $t_end)){
			#AND there's only a SINGLE liftover match in this region...
			my $nummatches= $regions{$cds3}{$intron3};
			unless ($nummatches == 1){
				next;
			}
			$checked++;						
			my $t_sflush= $target{$cds3}{$intron3}{'start_flush'};
			my $t_eflush= $target{$cds3}{$intron3}{'end_flush'};

			#######CHECK INTRON START POSITIONS#######
			#if both the query & stored target starts are flush...expect target start to match liftover start...
			#BUT- won't necessarily be true if there are multiple liftover matches all within same region...
			if (($q_sflush eq 'yes') && ($t_sflush eq 'yes')){
				#double-check liftover start coord matches stored expected target start
				unless ($i_start == $t_start){
					#if it doesn't match...update to the intron liftover start if it's a higher value
					#is more conservative, in case there is extra coding seq in target relative to galGal ref...
					if ($i_start > $t_start){
						$target{$cds3}{$intron3}{'start'}= $i_start;
						$updated_start++;
					}
				}
			}
			#if the query start is flush, but stored target start is not, update from liftover start...
			elsif (($q_sflush eq 'yes') && ($t_sflush eq 'no')){
				#overwrite the target intron start and flush stored info.
				$target{$cds3}{$intron3}{'start'}= $i_start;
				$target{$cds3}{$intron3}{'start_flush'}= 'yes';
				$updated_start++;
			}
			#if the query start and stored target start are not flush, update the start coord and leave 'flush' flag as is
			elsif (($q_sflush eq 'no') && ($t_sflush eq 'no')){
				$target{$cds3}{$intron3}{'start'}= $i_start;
				$updated_start++;
			}
			#if the query is not flush, but the target was...
			elsif (($q_sflush eq 'no') && ($t_sflush eq 'yes')){
				#just check that the stored start is actually smaller than liftover value...
				unless ($t_start <= $i_start){
					print "Error!  $cds3 $intron3 'expected' flush start: $t_start is not smaller than liftover non-flush start: $i_start\n";
					die;
				}
			}
			else{
				print "Unexpected 'start' comparison: $cds3 $intron3 query flush: $q_sflush target flush: $t_sflush\n";
				die;
			}

			#######CHECK INTRON END POSITIONS#######
			#if both the query & stored target ends are flush...expect target end to match liftover end...
			#BUT- won't necessarily be true if there are multiple liftover matches all within same region...
			if (($q_eflush eq 'yes') && ($t_eflush eq 'yes')){
				#double-check liftover start coord matches stored expected target start
				unless ($i_end == $t_end){
					#BUT- if it doesn't match...update to the intron liftover end if it's a lower value
					#is more conservative, in case there is extra coding seq in target relative to galGal ref...
					if ($i_end < $t_end){
						$target{$cds3}{$intron3}{'end'}= $i_end;
						$updated_end++;
					}
				}
			}
			#if the query end is flush, but stored target end is not, update from liftover end...
			elsif (($q_eflush eq 'yes') && ($t_eflush eq 'no')){
				#overwrite the target intron end and flush stored info.
				$target{$cds3}{$intron3}{'end'}= $i_end;
				$target{$cds3}{$intron3}{'end_flush'}= 'yes';
				$updated_end++;
			}
			#if the query end and stored target end are not flush, update the end coord and leave 'flush' flag as is
			elsif (($q_eflush eq 'no') && ($t_eflush eq 'no')){
				$target{$cds3}{$intron3}{'end'}= $i_end;
				$updated_end++;
			}
			#if the query is not flush, but the target was...
			elsif (($q_eflush eq 'no') && ($t_eflush eq 'yes')){
				#just check that the stored end is actually larger than liftover value...
				unless ($t_end >= $i_end){
					print "Error!  $cds3 $intron3 'expected' flush end: $t_end is not larger than liftover non-flush end: $i_end\n";
					die;
				}
			}
			else{
				print "Unexpected 'end' comparison: $cds3 $intron3 query flush: $q_eflush target flush: $t_eflush\n";
				die;
			}
		}
	}
}
close(LIFT2);
print "\nFinished checking: $checked introns\n";
print "Updated: $updated_start start coordinates\n";
print "Updated: $updated_end end coordinates\n\n";

#output to .bed file with updated coordinates
open (BED, ">$spp\_nonoverlapping_introns.bed") || die "Cannot open output bedfile for writing: $!\n";

#also...store the count of how many target matching scaffolds there are...want to only keep 'single' cases...
my $matchfile= 'halLiftover/' . $spp . '_coding_exons_scaffold_match_info';
open (MATCHES, "<$matchfile") || die "Cannot open $matchfile for reading: $!\n";
my %goodmatch;
while(<MATCHES>){
	chomp($_);
	if ($_=~ /^(CDS)/){
		next;
	}
	else{
		my @sp= split(/\t/, $_);
		if ($sp[2] == 1){
			$goodmatch{$sp[0]}= 'exists';
		}
	}
}
close(MATCHES);

my $output= 0;
my $output2= 0;
my $fil= 0;
my $filtered2= 0;
my $filtered3= 0;
for my $key1 (sort {$a cmp $b} (keys %target)){
	for my $key2 (sort {$a <=> $b} (keys %{$target{$key1}})){
		my $len= $target{$key1}{$key2}{'end'} - $target{$key1}{$key2}{'start'};
		#will only have galGal length stored if this is a 'nonoverlapping' intron...
		my $glen= 'n/a';
		if (exists $gal{$key1}{$key2}){
			$glen= $gal{$key1}{$key2}{'end'} - $gal{$key1}{$key2}{'start'};
		}
		my $filtered= 'no';
		if ($glen eq 'n/a'){
			$filtered= 'yes';
			$fil++;
		}
		#filter unless this is a 'good' CDS with only a single target liftover region from chicken...
		unless (exists $goodmatch{$key1}){
			$filtered= 'yes';
			$fil++;
		}
		#filter any target introns larger than 100kb
		if ($len >= 100000){
			$filtered= 'yes';
			$filtered3++;
		}
		#for any target introns larger than 10kb, filter if target is more than 50% longer than galGal ref
		if (($len >= 10000) and ($glen ne 'n/a')){
			my $diff= (($len - $glen)/$len)*100;
			if ($diff > 50){
				$filtered= 'yes';
				$filtered3++;
			}
		}
		#for any target introns smaller than 100bp, filter if target is less than 50% of galGal ref
		if (($len <= 100) and ($glen ne 'n/a')){
			if ($len <= 0){
				$filtered= 'yes';
				$filtered2++;
			}
			else{
				my $diff2= (($len - $glen)/$len)*100;
				if ($diff2 < 50){
					$filtered= 'yes';
					$filtered2++;
				}
			}
		}

		#only output this intron if it was a 'nonoverlapping' one in galGal (e.g. doesn't overlap any exonic feature)
		#and is not deemed overly (suspiciously) long, or too short
		if ((exists $gal{$key1}{$key2}) and ($filtered eq 'no')){		
			print BED "$target{$key1}{$key2}{'scaffold'}\t$target{$key1}{$key2}{'start'}\t$target{$key1}{$key2}{'end'}\tCDS=$key1,intron=$key2\t0\t$target{$key1}{$key2}{'strand'}\n";
			$output2++;
		}		
		$output++;
	}
}
close(BED);
print "Output coordinates for: $output target introns\n\n";
print "Filtered: $fil introns (not 'nonoverlapping', or with multiple target matches)\n";
print "Filtered: $filtered3 introns (too big)\n";
print "Filtered: $filtered2 introns (too small)\n";
print "Output .bed coordinates for: $output2 'NONOVERLAPPING' introns\n\n";
