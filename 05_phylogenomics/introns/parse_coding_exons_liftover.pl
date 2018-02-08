#!/usr/bin/perl
use strict;
use warnings;

#specify species on command line (e.g. droNov)
my $spp= shift(@ARGV);

#store info. for galGal reference coding exons
my $galfile= 'input_annotations/galGal4_coding_exons_info';
open(GAL, "<$galfile") || die "Cannot open $galfile for reading: $!\n";

my %gal_exons;
my %gal_ranks;
my $count= 0;
while(<GAL>){
	chomp($_);
	#skip header
	if ($_=~ /^(CDS)/){
		next;
	}
	else{
		my @split= split(/\t/, $_);
		my ($cds, $scaff, $start, $end, $strand, $rank)= ($split[0], $split[5], $split[6], $split[7], $split[8], $split[10]);
	
		$gal_exons{$cds}{$start}{'strand'}= $strand;
		$gal_exons{$cds}{$start}{'rank'}= $rank;
		$count++;

		#also store keyed by cds, then rank
		$gal_ranks{$cds}{$rank}{'scaffold'}= $scaff;
		$gal_ranks{$cds}{$rank}{'strand'}= $strand;
		$gal_ranks{$cds}{$rank}{'start'}= $start;
		$gal_ranks{$cds}{$rank}{'end'}= $end;
	}
}
close(GAL);
my $scalar_cds= scalar(keys %gal_exons);
print "\nFinished storing: $count exons for $scalar_cds galGal ref CDS\n";

#parse galGal-> target spp halLiftover output
my $targetfile= 'halLiftover/' . $spp . '_coding_exons_liftover.psl';
my %target_exons;
open(TARGET, "<$targetfile") || die "Cannot open $targetfile for reading: $!\n";

my $count2= 0;
while(<TARGET>){
	chomp($_);
	my @split2= split(/\t/, $_);
	my $id= $split2[0];
	#parse the galGal cds, start & end from ID
	my $cds2;
	my $start2;
	my $end2;
	if ($id=~ /^(CDS=)(.+?)(_)(\d+?)(-)(\d+?)(,)/){
		$cds2= $2;
		$start2= $4;
		$end2= $6;
	}
	else{
		print "Error parsing: $id\n";
		die;
	}

	#get the expected strand & exon rank from the stored galGal info
	my $strand2= $gal_exons{$cds2}{$start2}{'strand'};
	my $rank2= $gal_exons{$cds2}{$start2}{'rank'};

	#parse the query/target strands from liftover output & double-check query (galGal ref) matches expected
	my $strands= $split2[9];
	my $qstrand;
	my $tstrand;
	if ($strands=~ /^([+-])([+-])$/){
		$qstrand= $1;
		$tstrand= $2;
	}
	else{
		print "Error parsing query/target strands from: $strands\n";
		die;
	}
	unless ($qstrand eq $strand2){
		print "Error!  Query strand: $qstrand does not match expected $strand2\n";
		die;
	}

	#get the query 'matched region' start & end, and record if match extends right to exon edges
	my $qstart= $split2[12];
	my $qend= $split2[13];
	my $start_flush= 'no';
	my $end_flush= 'no';

	#if the matched region start equals the annotated exon start, set 'start flush' to yes...
	if ($qstart == $start2){
		$start_flush= 'yes';
	}
	#if matched region end matches annotated exon end...
	if ($qend == $end2){
		$end_flush= 'yes';
	}

	#get info. for the target scaffold, start, and end
	my $tscaff= $split2[14];
	my $tstart= $split2[16];
	my $tend= $split2[17];

	#store info in %target_exons, keyed by cds, then scaffold, exon rank, and exon start
	$target_exons{$cds2}{$tscaff}{$rank2}{$tstart}{'end'}= $tend;
	$target_exons{$cds2}{$tscaff}{$rank2}{$tstart}{'strand'}= $tstrand;
	$target_exons{$cds2}{$tscaff}{$rank2}{$tstart}{'start_flush'}= $start_flush;
	$target_exons{$cds2}{$tscaff}{$rank2}{$tstart}{'end_flush'}= $end_flush;
	$count2++;
}
close(TARGET);
my $scalar_cds2= scalar(keys %target_exons);
print "\nFinished storing: $count2 exons for $scalar_cds2 target CDS halLiftover annotations...\n";

#if there's a single scaffold with liftovers...keep; if there are multiple, keep scaffold with greatest number of liftover exons
#and, allow only a single match for each exon within a given scaffold
my %goodmatch;
my %stored1;
my %stored2;
my %stored3;
for my $cdskey (sort {$a cmp $b} (keys %target_exons)){
	my %excount;
	my %exs;
	my $numscaffs= scalar(keys %{$target_exons{$cdskey}});
	$stored1{$cdskey}{'num_scaffs'}= $numscaffs;

	for my $scaffkey (sort {$a cmp $b} (keys %{$target_exons{$cdskey}})){
		$stored2{$cdskey}{$scaffkey}{'exon_ranks'}= '';
		my $numexons= scalar(keys %{$target_exons{$cdskey}{$scaffkey}});
		for my $rankkey (sort {$a <=> $b} (keys %{$target_exons{$cdskey}{$scaffkey}})){
			my $numseqs= scalar(keys %{$target_exons{$cdskey}{$scaffkey}{$rankkey}});
			unless ($numseqs == 1){
				#if this exon has multiple matches in this scaffold...skip...so decrement count of 'good' exons...
				$numexons--;
			}
			#if there is a single match...store this scaffold/exon rank as 'good'
			else{
				$exs{$cdskey}{$scaffkey}{$rankkey}= 'exists';
				$stored2{$cdskey}{$scaffkey}{'exon_ranks'}= $stored2{$cdskey}{$scaffkey}{'exon_ranks'} . $rankkey . '|';
			}
		}
		#and, store the final number of good exon matches for this scaffold...
		$excount{$numexons}{$scaffkey}= 'exists';
		$stored2{$cdskey}{$scaffkey}{'num_exons'}= $numexons;
		$stored3{$cdskey}{$numexons}= 'exists';
	}
	#after checking all scaffolds for this cds...
	#if there's a single scaffold with good exon matches...use
	my $scal1= scalar(keys %{$exs{$cdskey}});
	if ($scal1 == 1){
		#store all the good exons to use in %goodmatch
		for my $usescaff (keys %{$exs{$cdskey}}){
			for my $userank (keys %{$exs{$cdskey}{$usescaff}}){
				$goodmatch{$cdskey}{$usescaff}{$userank}= 'exists';
			}
		}
	}
	#if there are multiple scaffolds with good exon matches, use scaffold with highest # of matches, or skip if there's a tie for the highest #
	else{
		my @counts= sort {$a <=> $b} (keys %excount);
		my $highest= pop(@counts);
		my $scal2= scalar(keys %{$excount{$highest}});
		#if there's a single scaffold with this highest number of 'good' exons...
		if ($scal2 == 1){
			for my $usescaff2 (keys %{$excount{$highest}}){
				for my $userank2 (keys %{$exs{$cdskey}{$usescaff2}}){
					$goodmatch{$cdskey}{$usescaff2}{$userank2}= 'exists';
				}
			}
		}
	}
}
my $scalgood= scalar(keys %goodmatch);
print "\nStored 'good' scaffold for: $scalgood CDS\n";

#now, loop through a 2nd time & keep 'good' liftover entries only if exons are all on same strand and in expected order...
my %goodexons;
my $goodcount= 0;
my %badstrand;
my %badorder;
for my $cdskey2 (sort {$a cmp $b} (keys %target_exons)){
	for my $scaffkey2 (sort {$a cmp $b} (keys %{$target_exons{$cdskey2}})){
		#only keep if it's 'good'
		unless (exists $goodmatch{$cdskey2}{$scaffkey2}){
			next;
		}
		#otherwise, check that all the matching exons are on the same strand & in the expected order...
		else{
			my @starts;
			my $storedstrand;
			#check against the number of 'good' exons ranks stored
			my $numranks= scalar(keys %{$goodmatch{$cdskey2}{$scaffkey2}});
			for my $rankkey2 (sort {$a <=> $b} (keys %{$target_exons{$cdskey2}{$scaffkey2}})){
				#skip if this wasn't stored as a 'good' exon
				unless (exists $goodmatch{$cdskey2}{$scaffkey2}{$rankkey2}){
					next;
				}
				for my $startkey2 (keys %{$target_exons{$cdskey2}{$scaffkey2}{$rankkey2}}){
					my $checkstrand= $target_exons{$cdskey2}{$scaffkey2}{$rankkey2}{$startkey2}{'strand'};
					if (defined($storedstrand)){
						if ($checkstrand eq $storedstrand){
							push(@starts, $startkey2);
						}
						else{
							$badstrand{$cdskey2}= 'exists';
						}
					}
					#if it's the first exon we're looping through...
					else{
						push(@starts, $startkey2);
						$storedstrand= $checkstrand;
					}
				}
			}
			#now, as long as we stored the same # of start positions are there are exon ranks, strand is good
			my $scalarstarts= scalar(@starts);
			unless ($scalarstarts == $numranks){
				next;
			}
			else{
				#...and as long as exon starts are good...store the exons from this CDS as 'good'
				#if we're on the + strand in target sequence, expect start positions to increase as we go through the exon ranks...
				my $storedstart;
				my $bad= 'no';
				if ($storedstrand eq '+'){
					#set initial 'storedstart' to -1
					$storedstart= -1;
					for my $s1 (@starts){
						if ($s1 > $storedstart){
							#this is as expected, so just update storedstart
							$storedstart= $s1;
							next;
						}
						else{
							#unexpected, so set 'bad' flag to 'yes'
							$bad= 'yes';
							$badorder{$cdskey2}= 'exists';
							next;
						}
					}
				}
				#if we're on the - strand, expect starts to get smaller...
				elsif ($storedstrand eq '-'){
					#set initial 'storedstart' to arbitrarily high number that will always be > than scaffold length
					$storedstart= 999999999999999999;
					for my $s2 (@starts){
						if ($s2 < $storedstart){
							$storedstart= $s2;
							next;
						}
						else{
							$bad= 'yes';
							$badorder{$cdskey2}= 'exists';
							next;
						}
					}
				}
				else{
					print "Error: unexpected strand $storedstrand\n";
					die;
				}
				#now, as long as 'bad' flag is still 'no', we'll keep this cds...
				if ($bad eq 'no'){
					#store this 'good' scaffold that we're using...
					$stored1{$cdskey2}{'final_scaff'}= $scaffkey2;
					for my $rankkey3 (keys %{$target_exons{$cdskey2}{$scaffkey2}}){
						#again, skip if not 'good'
						unless (exists $goodmatch{$cdskey2}{$scaffkey2}{$rankkey3}){
							next;
						}
						for my $startkey3 (keys %{$target_exons{$cdskey2}{$scaffkey2}{$rankkey3}}){
							$goodexons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'end'}= $target_exons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'end'};
							$goodexons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'strand'}= $target_exons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'strand'};
							$goodexons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'start_flush'}= $target_exons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'start_flush'};
							$goodexons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'end_flush'}= $target_exons{$cdskey2}{$scaffkey2}{$rankkey3}{$startkey3}{'end_flush'};
							$goodcount++;
						}
					}
				}
			}
		}
	}
}
my $scalargood= scalar(keys %goodexons);
print "\nRetained: $goodcount 'good' target exons from $scalargood CDS\n";
my $scalarbadstrand= scalar(keys %badstrand);
print "Omitted: $scalarbadstrand CDS with 'bad' strands (annotated to both + and -)\n";
my $scalarbadorder= scalar(keys %badorder);
print "Omitted: $scalarbadorder CDS with 'bad' exon order\n";

#print info. for these good exons/CDS
open(OUT, ">$spp\_good_coding_exons_info") || die "Cannot open outfile for writing: $!\n";
print OUT "CDS\tExon rank\tgalGal scaffold\tgalGal start\tgalGal end\tgalGal strand\t$spp scaffold\t$spp start\t$spp end\t$spp strand\t$spp start flush\t$spp end flush\n";

for my $k1 (sort {$a cmp $b} (keys %goodexons)){
	for my $k2 (sort {$a cmp $b} (keys %{$goodexons{$k1}})){
		for my $k3 (sort {$a <=> $b} (keys %{$goodexons{$k1}{$k2}})){
			#get the galGal ref info. for this cds/exon rank
			my $gscaff= $gal_ranks{$k1}{$k3}{'scaffold'};
			my $gstart= $gal_ranks{$k1}{$k3}{'start'};
			my $gend= $gal_ranks{$k1}{$k3}{'end'};
			my $gstrand= $gal_ranks{$k1}{$k3}{'strand'};

			for my $k4 (keys %{$goodexons{$k1}{$k2}{$k3}}){
				my $goodend= $goodexons{$k1}{$k2}{$k3}{$k4}{'end'};
				my $goodstrand= $goodexons{$k1}{$k2}{$k3}{$k4}{'strand'};
				my $goodsflush= $goodexons{$k1}{$k2}{$k3}{$k4}{'start_flush'};
				my $goodeflush= $goodexons{$k1}{$k2}{$k3}{$k4}{'end_flush'};

				print OUT "$k1\t$k3\t$gscaff\t$gstart\t$gend\t$gstrand\t$k2\t$k4\t$goodend\t$goodstrand\t$goodsflush\t$goodeflush\n";
			}
		}
	}
}
close(OUT);

#and...print out some info. about how many possible scaffolds had matches, which was used, etc.
open(OUT2, ">$spp\_coding_exons_scaffold_match_info") || die "Cannot open outfile for writing: $!\n";
print OUT2 "CDS\t#galGal exons\t# Scaffolds with matches\tFinal scaffold used\tHighest # exon matches\tNext highest\tScaffold info. (ID, # 'good' exons, exon ranks for each candidate scaffold)\n";

for my $ka (sort {$a cmp $b} (keys %stored1)){
	my $string= $stored1{$ka}{'final_scaff'};
	unless(defined($string)){
		$string= 'n/a';
	}
	my $galnum= scalar(keys %{$gal_ranks{$ka}});
	print OUT2 "$ka\t$galnum\t$stored1{$ka}{'num_scaffs'}\t$string\t";
	my @nums= sort {$a <=> $b} (keys %{$stored3{$ka}});
	my $highest= pop(@nums);
	unless(defined($highest)){
		$highest= 'n/a';
	}
	my $nexthigh= pop(@nums);
	unless(defined($nexthigh)){
		$nexthigh= 'n/a';
	}
	print OUT2 "$highest\t$nexthigh";
	for my $kb (sort {$a cmp $b} (keys %{$stored2{$ka}})){
		
		print OUT2 "\t$kb\t$stored2{$ka}{$kb}{'num_exons'}\t$stored2{$ka}{$kb}{'exon_ranks'}";
	}
	print OUT2 "\n";
}
close(OUT2);		
