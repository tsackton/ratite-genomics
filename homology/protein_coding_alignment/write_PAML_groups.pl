#!/usr/bin/perl
use strict;
use warnings;

#list of all taxa (outgroups, neognaths, palaeognaths)
my @taxa= qw( allMis anaPla anoCar aptFor aptHaa aptOwe aptRow aquChr balReg calAnn casCas chaPel chaVoc chrPic colLiv corBra cryCin cucCan droNov egrGar eudEle falPer ficAlb fulGla galGal geoFor halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen serCan strCam taeGut tinGut );

#list of all palaeognath taxa
my @palaeo_list= qw( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );
#map array entries to hash
my %palaeos = map { $_ => 1 } @palaeo_list;

#path to original (unfiltered) alignment summary stats
my $summaryfile= '/home/alison/Desktop/oma_duplicate_genes/alignment_stats/fastas_aligned/alignment_summary_stats';
open(SUMMARY1, "<$summaryfile") || die "Cannot open $summaryfile for reading: $!\n";
#store the average sequence length for each locus in %avglen
my %avglen;

while(<SUMMARY1>){
	chomp($_);
	#skip header 
	if ($_=~ /^(Batch)/){
		next;
	}
	my @split1= split(/\t/, $_);
	my ($hog1, $len1)= ($split1[1], $split1[4]);
	$avglen{$hog1}= $len1;
}
close(SUMMARY1);
my $scalar1= scalar(keys %avglen);
print "\nFinished storing average aligned sequence length for: $scalar1 loci\n";

#now, for each taxon (except outgroups), ID any seqs to omit from final group
#exclusion criteria: if length of filtered seq is less than 50% of the original input length, 
#OR if the unfiltered length is less than 50% of the average input sequence length of all taxa (in original unfiltered alignment)
#OR if filtered seq has > 1gap/bp aligned sequence
my %remove;
my %goodseqs;
my %allseqs;
my %badseqs;
my %dontcount;
for my $spp (@taxa){
	my $count= 0;
	if (($spp eq 'allMis') || ($spp eq 'anoCar') || ($spp eq 'chrPic')){
		next;
	}
	else{
		print "\nProcessing alignment stats for: $spp\n";
		#path to the original (unfiltered) alignment stats
		my $orig= '/home/alison/Desktop/oma_duplicate_genes/alignment_stats/fastas_aligned/' . $spp . '_alignment_stats_edited';
		#path to the filtered alignment stats
		my $filter= '/home/alison/Desktop/oma_duplicate_genes/alignment_stats/fastas_aligned_filter3/' . $spp . '_alignment_stats_edited';

		#store all filtered lengths in temp %flen hash for this species, & gaps/bp in %fgaps 
		my %flen;
		my %fgaps;
		open(FILTER, "<$filter") || die "Cannot open $filter for reading: $!\n";
		while(<FILTER>){
			chomp($_);
			if ($_=~ /^(Batch)/){
				next;
			}
			else{
				my @split2= split(/\t/, $_);
				#skip if no seqs. for this spp. in this locus
				if ($split2[2] eq 'n/a'){
					next;
				}
				else{
					#store filtered length, keyed by locus then seqid (to avoid overwriting if multiple spp. seqs in same group)
					$flen{$split2[1]}{$split2[2]}= $split2[3];
					#store gaps per bp aligned
					$fgaps{$split2[1]}{$split2[2]}= $split2[4];
				}
			}
		}
		close(FILTER);

		#now, loop through unfiltered stats...
		open(ORIG, "<$orig") || die "Cannot open $orig for reading: $!\n";
		while(<ORIG>){
			chomp($_);
			if ($_=~ /^(Batch)/){
				next;
			}
			else{
				my @split3= split(/\t/, $_);
				#skip if no seqs. for this spp. in this locus
				if ($split3[2] eq 'n/a'){
					next;
				}
				else{
					my $hog3= $split3[1];
					my $id3= $split3[2];
					#original length
					my $origlen= $split3[3];
					#get the average input length for all spp for this locus...
					my $getavg= $avglen{$hog3};
					my $checkavg= ($origlen/$getavg)*100;

					#get the filtered seq length for this sequence
					my $getlen;
					if (exists $flen{$hog3}{$id3}){
						$getlen= $flen{$hog3}{$id3};
					}
					#if sequence was entirely filtered out (won't have appeared in filter3 file...)
					else{
						$getlen= 0;
						#but make sure we don't count this in calculating # of good spp below...
						$dontcount{$spp}{$hog3}{$id3}= 'exists';
					}
					my $checklen= ($getlen/$origlen)*100;

					#if orig seq. is less than 50% of avg. unfiltered seq length OR filtered length is less than 50% of orig...we'll remove seq
					if (($checkavg < 50) || ($checklen < 50)){
						$remove{$spp}{$hog3}{$id3}= 'exists';
						$badseqs{$hog3}{$spp}{$id3}= 'exists';
					}
					#also, remove seq if we have > 1gap/bp aligned in filtered dataset
					if (exists $fgaps{$hog3}{$id3}){
						if ($fgaps{$hog3}{$id3} > 1){
							$remove{$spp}{$hog3}{$id3}= 'exists';
							$badseqs{$hog3}{$spp}{$id3}= 'exists';
						}
					}
					#and, store ALL seqs...
					$allseqs{$hog3}{$spp}{$id3}= 'exists';
				}
			}
		}
		close(ORIG);
		my $scalar_remove= scalar(keys %{$remove{$spp}});
		print "Remove: $scalar_remove $spp seqs (too short/too much filtered/too gappy)\n";
	}
}

#now, loop through filtered alignment summary stats
#keep loci if: no more than 50% of all bird spp. missing AND no more than 50% of all paleognath spp. missing
#then, as a reasonable threshold for allowable duplications: no more than 3 'good' seqs/species and tot # seqs not more than 1.5X tot species in alignment
my %goodhogs;

#path to filtered alignment summary stats
my $summaryfile2= '/home/alison/Desktop/oma_duplicate_genes/alignment_stats/fastas_aligned_filter3/alignment_summary_stats';
open(SUMMARY2, "<$summaryfile2") || die "Cannot open $summaryfile2 for reading: $!\n";
my $tot= 0;
my $goodnum= 0;
print "\nProcessing: HOGs...\n";
while(<SUMMARY2>){
	chomp($_);
	#skip header 
	if ($_=~ /^(Batch)/){
		next;
	}
	else{
		$tot++;
		if ($tot=~ /^(\d+)(0{3})$/){
			print "Processed: $tot\n";
		}
		my %hoginfo;
		my @hogvals;
		my $sppcount= 0;
		my $palaeocount= 0;
		my @spl= split(/\t/, $_);
		my $grp= $spl[1];
		my @tax2= @taxa;
		my $totgoodseqs= 0;

		#store the number of seqs for each spp.
		for (my $i= 8; $i <= 49; $i++){
			#get the spp. name for this entry
			my $taxid= shift(@tax2);
			#skip outgroups
			if (($taxid eq 'allMis') || ($taxid eq 'anoCar') || ($taxid eq 'chrPic')){
				next;
			}
			my $numseqs= $spl[$i];
			#replace n/a with 0
			if ($numseqs eq 'n/a'){
				$numseqs= 0;
			}
			#get number of bad seqs. to remove for this spp.
			my $todelete= scalar(keys %{$remove{$taxid}{$grp}});
			#BUT- don't count any that were entirely filtered out b/c won't have appeared in filter3 tally to begin with...
			my $toadd= scalar(keys %{$dontcount{$taxid}{$grp}});
#			print "$taxid: delete $todelete\n";
			#add entirely filtered out seqs ONLY if tally was nonzero to begin with...
			if ($numseqs != 0){
				$numseqs= $numseqs + $toadd;
			}
			#then, delete all bad seqs.
			$numseqs= $numseqs - $todelete;

			#if we have at least 1 good seq for this spp...	
			if ($numseqs >= 1){
				$sppcount++;
				#if it's a palaeognath
				if (exists $palaeos{$taxid}){
					$palaeocount++;
				}
				#and, store the number of good seqs for this spp, keyed by number of seqs then by taxa that have that number...
				$hoginfo{$numseqs}{$taxid}= 'exists';
				#and, push value into @hogvals
				push(@hogvals, $numseqs);
				$totgoodseqs+= $numseqs;
			}
		}
		#now that we have the number of good seqs for each (non-outgroup) species...
		#keep HOG if: we have at least 50% of birds, and at least 50% of palaeos
		if (($sppcount > 19) && ($palaeocount >= 6)){
			$goodnum++;
			#criteria for allowable gene duplication criteria...
			my @sorted_vals= sort {$a <=> $b} (@hogvals);
			my $lowest= $sorted_vals[0];
			my $lastindex= $#sorted_vals;
			my $highest= $sorted_vals[$lastindex];

			#Tim's suggestion: max 3 seqs/spp, no restriction on # spp with duplications, but tot #seqs not more than 1.5X the number of spp
			my $ratio= $totgoodseqs/$sppcount;
			if (($highest <= 3) && ($ratio <= 1.5)){
				$goodhogs{$grp}= $totgoodseqs;
			}
		}
	}
}
close(SUMMARY2);
my $scalargood= scalar(keys %goodhogs);
#print "\nKeep: $scalargood 'good' HOGs (of $tot total)\n";
print "\nTotal HOGs: $tot\n";
print "HOGs meeting min. species/seq. length requirements: $goodnum\n";
print "HOGs also meeting max. duplication cutoff: $scalargood\n\n";

#so, for these good groups, output:
#the list of 'good' HOGs to use, with the seqs to keep for each HOG (taxIDs & protein IDs)
open(OUT, ">good_PAML_HOG_protein_info") || die "Cannot open outfile for writing: $!\n";
print OUT "HOG\tSpecies\tProtein ID\n";
my $output= 0;
for my $key1 (sort {$a cmp $b} (keys %goodhogs)){
	$output++;
	#double-check that we output the expected # of 'good' seqs for each group...
	my $checkseqs= $goodhogs{$key1};
	my $outputseqs= 0;
	for my $key2 (sort {$a cmp $b} (keys %{$allseqs{$key1}})){
		for my $key3 (sort {$a cmp $b} (keys %{$allseqs{$key1}{$key2}})){
			#only print if not stored as a bad seq...
			unless(exists $badseqs{$key1}{$key2}{$key3}){
				#parse just the actual protein ID
				my $protid;
				if ($key3=~ /^($key2)(_)(.+)$/){
					$protid= $3;
				}
				else{
					print "Error parsing protein ID from: $key3\n";
					die;
				}
				print OUT "$key1\t$key2\t$protid\n";
				$outputseqs++;
			}
		}
	}
	unless ($outputseqs == $checkseqs){
		print "Error!  Output $outputseqs good seqs for $key1 (expected $checkseqs)\n";
		print "All seqs:\n";
		for my $k2 (sort {$a cmp $b} (keys %{$allseqs{$key1}})){
			for my $k3 (sort {$a cmp $b} (keys %{$allseqs{$key1}{$k2}})){
				print "$k3\n";
			}
		}
		print "Bad seqs:\n";
		for my $k4 (sort {$a cmp $b} (keys %{$badseqs{$key1}})){
			for my $k5 (sort {$a cmp $b} (keys %{$badseqs{$key1}{$k4}})){
				print "$k5\n";
			}
		}
		die;
	}
}
print "Ouput protein ID info for: $output HOGs\n\n";
