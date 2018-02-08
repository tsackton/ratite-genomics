#!/usr/bin/perl
use strict;
use warnings;

#get input species from command line (e.g. droNov)
my $spp= shift(@ARGV);

#store info from parsing of coding exons liftover
my $exoninfo= $spp . '_good_coding_exons_info';
open(EXINFO, "<$exoninfo") || die "Cannot open $exoninfo for reading: $!\n";

my %exons;
while(<EXINFO>){
	chomp($_);
	if ($_=~ /^(CDS)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		my $cds= $split[0];
		my $rank= $split[1];
		
		$exons{$cds}{$rank}{'scaffold'}= $split[6];
		$exons{$cds}{$rank}{'start'}= $split[7];
		$exons{$cds}{$rank}{'end'}= $split[8];
		$exons{$cds}{$rank}{'strand'}= $split[9];
		$exons{$cds}{$rank}{'start_flush'}= $split[10];
		$exons{$cds}{$rank}{'end_flush'}= $split[11];
	}
}
close(EXINFO);

open(OUT, ">$spp\_expected_intron_coords") || die "Cannot open outfile for writing: $!\n";
print OUT "CDS\tIntron\tScaffold\tStrand\tStart\tEnd\n";
for my $cdskey (sort {$a cmp $b} (keys %exons)){
	for my $rankkey (sort {$a <=> $b} (keys %{$exons{$cdskey}})){
		#if we also have the next exon rank stored...
		my $rankkey2= $rankkey + 1;
		if (exists $exons{$cdskey}{$rankkey2}){
			my $intnum= $rankkey;
			
			my $scaffold= $exons{$cdskey}{$rankkey}{'scaffold'};
			my $strand= $exons{$cdskey}{$rankkey}{'strand'};

			my $start1= $exons{$cdskey}{$rankkey}{'start'};
			my $end1= $exons{$cdskey}{$rankkey}{'end'};
			my $sflush1= $exons{$cdskey}{$rankkey}{'start_flush'};
			my $eflush1= $exons{$cdskey}{$rankkey}{'end_flush'};

			my $start2= $exons{$cdskey}{$rankkey2}{'start'};
			my $end2= $exons{$cdskey}{$rankkey2}{'end'};
			my $sflush2= $exons{$cdskey}{$rankkey2}{'start_flush'};
			my $eflush2= $exons{$cdskey}{$rankkey2}{'end_flush'};

			my $intron_start;
			my $intron_end;

			#if this cds is on the + strand...
			if ($strand eq '+'){
				#if the lower ranked exon has a flush end, the intron start is that end (e.g. that exon end [1-based coord] is the intron start in 0-based coords)
				if ($eflush1 eq 'yes'){
					$intron_start= $end1;
				}
				#if it's not flush, the intron start is greater than that end...
				else{
					$intron_start= '>' . $end1;
				}
				#if the higher ranked exon has a flush start, the intron end (in 1-based coords) is that start (0-based)
				if ($sflush2 eq 'yes'){
					$intron_end= $start2;
				}
				#if not flush, the intron end is less than that start...
				else{
					 $intron_end= '<' . $start2;
				}
			}
			#if this cds is on the - strand...
			elsif ($strand eq '-'){
				if ($sflush1 eq 'yes'){
					#NB- 'start' position must always be smaller than 'end', so call this the 'intron_end'
					$intron_end= $start1;
				}
				else{
					$intron_end= '<' . $start1;
				}
				if ($eflush2 eq 'yes'){
					$intron_start= $end2;
				}
				else{
					$intron_start= '>' . $end2;
				}
			}
			else{
				print "CDS: $cdskey, Unexpected strand: $strand\n";
				die;
			}
			print OUT "$cdskey\t$intnum\t$scaffold\t$strand\t$intron_start\t$intron_end\n";
		}
	}
}
close(OUT);			
