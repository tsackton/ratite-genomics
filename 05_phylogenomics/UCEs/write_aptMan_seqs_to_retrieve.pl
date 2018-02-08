#!/usr/bin/perl
use strict;
use warnings;

use Cwd;

my $queryspp= 'aptMan';

my @hitspp= qw ( aptHaa aptOwe aptRow );


#store the 'good' blastn hits from each query species...
my $count= 0;
my $kept= 0;
my %hits;
for my $taxon (@hitspp){
	#my $infile= "$queryspp-$taxon\_concatenated_blastn_results";
	my $infile= "$queryspp-$taxon\_blastn_results_all";
	print "\nProcessing: $infile\n";

	open(INFILE, "<" . $infile) || die "Cannot open $infile for reading: $!\n";
	while(<INFILE>){
		chomp($_);
		#if ($_ =~ /^(Batch)/){
		if ($_ =~ /^(Query)/){
			next;
		}
		else{
			$count++;
			my @split= split(/\t/, $_);
			#skip if no hits
			if ($split[2] eq 'NO HITS'){
				next;
			}

			#my $batch= $split[0];
			#my $query= $split[1];
			#my $qlen= $split[2];
			#my $hitname= $split[4];
			#my $hitstrand= $split[7];
			#my $eval= $split[8];
			#my $hsps= $split[9];
			#my $aln_qlen= $split[10];
			#my $fracident= $split[13];
			#my $hitstart= $split[17];
			#my $hitend= $split[18];

			#if not done with 'batches'...
			my $query= $split[0];
			my $qlen= $split[1];
			my $hitname= $split[3];
			my $hitstrand= $split[6];
			my $eval= $split[7];
			my $hsps= $split[8];
			my $aln_qlen= $split[9];
			my $fracident= $split[12];
			my $hitstart= $split[16];
			my $hitend= $split[17];

			#only keep as a 'good' hit if: single HSP, evalue <= 1e-10, covers 50% of query, frac identical >= 0.8
			my $coverage= ($aln_qlen/$qlen)*100;

			if (($hsps == 1) and ($eval <= 1e-10) and ($coverage >= 50) and ($fracident >= 0.8)){
				$kept++;
				#store in %hits, keyed by query name, then hit name, then strand & with batch, query length hit start and hit end as values
				#$hits{$query}{$hitname}{$hitstrand}{'batch'}= $batch;
				$hits{$query}{$hitname}{$hitstrand}{'query_length'}= $qlen;

				#keep the lowest hit start & highest hit end
				if (exists $hits{$query}{$hitname}{$hitstrand}{'start'}){
					if ($hitstart < $hits{$query}{$hitname}{$hitstrand}{'start'}){
						$hits{$query}{$hitname}{$hitstrand}{'start'}= $hitstart;
					}
				}
				else{
					$hits{$query}{$hitname}{$hitstrand}{'start'}= $hitstart;
				}
				if (exists $hits{$query}{$hitname}{$hitstrand}{'end'}){
					if ($hitend > $hits{$query}{$hitname}{$hitstrand}{'end'}){
						$hits{$query}{$hitname}{$hitstrand}{'end'}= $hitend;
					}
				}
				else{
					$hits{$query}{$hitname}{$hitstrand}{'end'}= $hitend;
				}
			}
		}
	}
	close(INFILE);
	print "Finished processing: $count hits for $taxon\n";
	print "Kept: $kept 'good' hits\n\n";
	$count= 0;
	$kept= 0;
}

#now, keep 'really' good hits if we've stored a single hit scaffold on the same strand for a given query, and the start & end coordinates aren't too much bigger than query
my $outfile= "$queryspp\_seqs_to_retrieve";
open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
print OUT "Query\tQuery length\tHit\tHit strand\tHit start\tHit end\n";

my $qcount= 0;
my $goodcount= 0;
for my $querykey (sort {$a cmp $b} (keys %hits)){
	$qcount++;
	my $scalarhits= scalar(keys %{$hits{$querykey}});
	if ($scalarhits == 1){
		for my $hitkey (keys %{$hits{$querykey}}){
			my $scalarstrands= scalar(keys %{$hits{$querykey}{$hitkey}});
			if ($scalarstrands == 1){
				for my $strandkey (keys %{$hits{$querykey}{$hitkey}}){
					my $qlen_val= $hits{$querykey}{$hitkey}{$strandkey}{'query_length'};
					#my $batch_val= $hits{$querykey}{$hitkey}{$strandkey}{'batch'};
					
					my $start_val= $hits{$querykey}{$hitkey}{$strandkey}{'start'};
					my $end_val= $hits{$querykey}{$hitkey}{$strandkey}{'end'};

					my $hitlen= ($end_val - $start_val) + 1;
					#if hit length is more than 10% longer than the input query length...omit
					my $checklen= (($hitlen - $qlen_val)/$qlen_val)*100;
					if ($checklen <= 110){
						print OUT "$querykey\t$qlen_val\t$hitkey\t$strandkey\t$start_val\t$end_val\n";
						$goodcount++;
					}
				}
			}
		}
	}
}
close(OUT);
print "Finished processing: $qcount queries with 'good' hit(s)\n";
print "Retained: $goodcount seqs. to retrieve\n\n";
