#!/usr/bin/perl
use strict;
use warnings;

#pass target species code via command line (e.g. aptOwe)
my $spp= shift(@ARGV);

#reference species used for halLiftover
my @refspp= qw ( galGal strCam tinGut );

#want to have final alignments all being on + strand in galGal reference, so store...
my $galbed= 'galGal_uces.bed';
my %galstrand;
open (GALBED, "<$galbed") || die "Cannot open $galbed for reading: $!\n";
while(<GALBED>){
	chomp($_);
	my @spl= split(/\t/, $_);
	if ($spl[3]=~ /^(galGal,)(.+)$/){
		$galstrand{$2}= $spl[5];
	}
	else{
		print "Error parsing galGal id: $spl[3]\n";
		die;
	}
}
close(GALBED);

my %loci;
for my $ref (@refspp){
	my $index= 0;
	print "\nParsing: $spp-$ref liftover...\n";
	my $liftfile= "halLiftover/$ref-$spp\_uces_halLiftover.psl";
	unless (-e $liftfile){
		if (($spp eq 'strCam') or ($spp eq 'tinGut') or ($spp eq 'galGal')){
			print "n/a...skipping!\n";
			next;
		}
		print "Error!!! No file: $liftfile\n";
		die;
	}
	open(LIFT, "<$liftfile") || die "Cannot open $liftfile for reading: $!\n";
	while(<LIFT>){
		$index++;
		chomp($_);
		my @split= split(/\t/, $_);
		my ($id, $scaff, $start, $end, $strand) = ($split[0], $split[14], $split[16], $split[17], $split[9]);
		my $locus;
		if ($id=~ /^([A-Za-z]{6})(,)(.+)$/){
			$locus= $3;
		}
		else{
			print "Error parsing locus from: $id\n";
			die;
		}
		$loci{$locus}{$ref}{$index}{'scaffold'}= $scaff;
		$loci{$locus}{$ref}{$index}{'start'}= $start;
		$loci{$locus}{$ref}{$index}{'end'}= $end;
		if ($strand=~ /^([+-])([+-])$/){
			$loci{$locus}{$ref}{$index}{'strand'}= $2;
		}
		else{
			print "Error parsing strand from: $strand\n";
			die;
		}
	}
	close(LIFT);
}

my %keep;
for my $key1 (sort {$a cmp $b} (keys %loci)){
	#only want to keep as a 'good' liftover if there is a unique liftover region for each reference (but allow locus to be missing from some reference taxa liftovers)
	my $flag= 'yes';
	for my $ck (keys %{$loci{$key1}}){
		my $ckscalar= scalar(keys %{$loci{$key1}{$ck}});
		if ($ckscalar > 1){
			$flag= 'no';
		}
	}
	#if we only had a single liftover region for each reference...
	if ($flag eq 'yes'){
		#check that all liftovers are on same scaffold & strand, and if so, use smallest start & largest end as long as they're all within a reasonable distance (500 bp)
		my $scaff;
		my $strand;
		my $start;
		my $end;
		my $good= 'yes';
		for my $key2 (keys %{$loci{$key1}}){
			for my $key3 (keys %{$loci{$key1}{$key2}}){
				my $getscaff= $loci{$key1}{$key2}{$key3}{'scaffold'};
				my $getstrand= $loci{$key1}{$key2}{$key3}{'strand'};
				my $getstart= $loci{$key1}{$key2}{$key3}{'start'};
				my $getend= $loci{$key1}{$key2}{$key3}{'end'};
				#if this is first loop, just store values...
				unless(defined($scaff)){
					$scaff= $getscaff;
					$strand= $getstrand;
					$start= $getstart;
					$end= $getend;
				}
				else{
					unless (($getscaff eq $scaff) && ($getstrand eq $strand)){
						$good= 'no';
					}
					else{
						my $startdiff= $start - $getstart;
						$startdiff= abs($startdiff);
						my $enddiff= $getend - $end;
						$enddiff= abs($enddiff);
						if (($getstart < $start) && ($startdiff <= 500)){
							$start= $getstart;
						}
						if (($getend > $end) && ($enddiff <= 500)){
							$end= $getend;
						}
					}
				}
			}
		}
		#if 'good' flag is still set to good, store the 'final' coords for this locus
		if ($good eq 'yes'){
			$keep{$key1}{'scaffold'}= $scaff;
			#check against the galGal annotation strand
			my $galstr= $galstrand{$key1};
			my $usestrand;
			if (defined($galstr)){
				if ($galstr eq '+'){
					$usestrand= $strand;
				}
				elsif ($galstr eq '-'){
					if ($strand eq '+'){
						$usestrand= '-';
					}
					elsif ($strand eq '-'){
						$usestrand= '+';
					}
					else{
						print "Locus: $key1\n";
						print "Error!  Unexpected strand: $strand\n";
						die;
					}
				}
				else{
					print "Error!  Unexpected strand: $strand\n";
				}
			}
			else{
				$usestrand= 'n/a';
				print "Warning: unknown strand orientation relative to chicken for: $key1 (output as 'n/a'...fix!!!)\n";
			}
			$keep{$key1}{'strand'}= $usestrand;
			#minus 1 from start to make 0-based
			$start= $start - 1;
			$keep{$key1}{'start'}= $start;
			$keep{$key1}{'end'}= $end;
		}
	}
}
my $scal1= scalar(keys %loci);
my $scal2= scalar(keys %keep);
print "Finished parsing: $scal1 loci for $spp...keep $scal2\n";	
#output to .bed file for this species...
my $outfile= $spp . '_uces_parsed_liftover.bed';
open(OUTBED, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
my $output= 0;
my $omit= 0;
for my $ka (sort {$a cmp $b} (keys %keep)){
	#skip any with unknown strand orientation...(if lots, could look up manually...)
	if ($keep{$ka}{'strand'} eq 'n/a'){
		$omit++;
		next;
	}
	$output++;
	print OUTBED "$keep{$ka}{'scaffold'}\t$keep{$ka}{'start'}\t$keep{$ka}{'end'}\tID=$ka\t0\t$keep{$ka}{'strand'}\n";
}
close(OUTBED);
print "Omitted: $omit loci with unknown galGal reference strand\n";
print "Output .bed coordinates for: $output loci...\n\n";		
