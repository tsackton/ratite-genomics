#!/usr/bin/perl
use strict;
use warnings;

my $statsfile= 'alignment_summary_stats';
open(STATS, "<$statsfile") || die "Cannot open $statsfile for reading: $!\n";

my %stats;
while(<STATS>){
	chomp($_);
	if ($_=~ /^(Batch)/){
		next;
	}
	else{
		my @spl= split(/\t/, $_);
		my ($batch, $locus, $numseqs, $avglen, $ident, $gaps)= ($spl[0], $spl[1], $spl[3], $spl[4], $spl[5], $spl[6]);
		
		#only keep if: no more than 1 missing spp., avg ident > 70%, gaps/bp aligned total < 0.5
		if (($numseqs >= 14) && ($ident >= 70) && ($gaps < 0.5)){
			#parse locus & intron
			my $cds;
			my $intron;
			if ($locus=~ /^(cds)(\d+?)(_)(intron)(\d+)$/){
				$cds= $1 . $2;
				$intron= $4 . $5;
			}
			else{
				print "Error parsing cds & intron from: $locus\n";
				die;
			}
			$stats{$cds}{$avglen}{$intron}= $batch;
		}
	}
}
close(STATS);

open(OUT, ">introns_dataset_list") || die "Cannot open outfile for writing: $!\n";
print OUT "Batch\tLocus\n";

#keep the intron with the highest number of seqs, and greatest average input sequence length for that # seqs
for my $cdskey (sort {$a cmp $b} (keys %stats)){
	my @nums= sort {$a <=> $b} (keys %{$stats{$cdskey}});
	my $highest= pop(@nums);

	for my $intronkey (keys %{$stats{$cdskey}{$highest}}){
		my $val= $stats{$cdskey}{$highest}{$intronkey};
		print OUT "$val\t$cdskey\_$intronkey\n";
	}
}
close(OUT);	
