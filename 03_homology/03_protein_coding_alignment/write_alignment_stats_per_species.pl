#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Math::Round;
use Statistics::Basic qw(:all);
use Cwd;

#list of all taxa (incl. outgroups, neognaths, palaeognaths)
my @taxa= qw( allMis anaPla anoCar aptFor aptHaa aptOwe aptRow aquChr balReg calAnn casCas chaPel chaVoc chrPic colLiv corBra cryCin cucCan droNov egrGar eudEle falPer ficAlb fulGla galGal geoFor halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen serCan strCam taeGut tinGut );

#glob all of the alignments in a specified directory
#&&& pass the directory via commandline...
my $indir= shift(@ARGV);
chomp($indir);
my $basedir= dirname($indir);

#pass the species as a 2nd command line argument
my $spp= shift(@ARGV);

#check if any results already output...and then open outfile for writing/appending
my %donejobs;
my $outfile= $spp . '_alignment_stats';
if (-e $outfile){
	open(DONE, "<$outfile") || die "Cannot open $outfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(Batch)/){
			next;
		}
		else{
			my @spl= split(/\t/, $_);
			my $scalarsplit= scalar(@spl);
			if (($scalarsplit == 49) or ($spl[2] eq 'n/a')){
				#store 'done' seqids as well as loci (b/c some OMAs have multiple seqs for a given spp & want all to run...)
				$donejobs{$spl[1]}{$spl[2]}= 'exists';
			}
		}
	}
	close(DONE);
	open(OUT, ">>$outfile") || die "Cannot open $outfile for appending: $!\n";
}
else{
	open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
	print OUT "Batch\tLocus\tID\tLength (bp; no gaps)\tGaps per bp aligned\tAvg. % ident (all spp.)\tAvg. % ident (non-self spp.)";
	for my $spp2 (@taxa){
		 print OUT "\t$spp-$spp2 % ident";
	}
	print OUT "\n";
}

#now, go through each alignment (in each 'batch')
my $count= 0;
for (my $i= 1; $i <= 4; $i++){
	#glob all of the fastas in that batch
	my @glob= <$indir/batch$i/*.fa>;
	#print "Globbing: $indir\batch$i/*.fa\n";
	
	my $scalarglob= scalar(@glob);
	print "Processing: $scalarglob fastas in batch$i for $spp\n";
	for my $filepath (sort {$a cmp $b} (@glob)){
		$count++;
		my $basename= basename($filepath);
		my $locus;
		if ($basename=~ /^(HOG2)(_)(\d+?)(_)/){
			$locus= $1 . $2 . $3;
		}
		else{
			print "Error parsing locus from: $basename\n";
			die;
		}
		my %statvals;

		#go through the alignment...
		my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
		my %seqnames;
		my %sppnames;
		while (my $aln= $alignio_in-> next_aln()){
			if ($aln-> is_flush()){
				#store all of the seq IDs, keyed first by taxon name...
				foreach my $seqobj ($aln-> each_seq()){
					my $id= $seqobj-> display_id();
					#skip if already done
					if (exists $donejobs{$locus}{$id}){
						next;
					}
					if ($id=~ /^([A-Za-z]{6})(.+)$/){
						my $shortid= $1;
						$seqnames{$shortid}{$id}= 'exists';
						#if this is a sequence for the species 'of interest'...
						if ($shortid eq $spp){
							$sppnames{$id}= 'exists';
						}
					}
					else{
						print "Error parsing taxon id from: $id\n";
						die;
					}
				}
			}
			else{
				print "Error!  $locus: alignment is not flush!\n";
				die;
			}
			#if there are no 'query' species seqs.
			my $scalarq= scalar(keys %sppnames);
			if ($scalarq == 0){
				print OUT "batch$i\t$locus\tn/a\n";
				next;
			}
			#now, for each of the target species sequences...
			foreach my $sppkey (sort {$a cmp $b} (keys %sppnames)){
				print "Processing: $sppkey for $locus (count: $count)\n";
				my $seqobj1= $aln-> get_seq_by_id($sppkey);
				my $seq1= $seqobj1-> seq();
				#replace leading/trailing gaps
				if ($seq1=~ /^(-+?)([^-])(.+)$/){
					$seq1= $2 . $3;
				}
				if ($seq1=~ /^(.+)([^-])(-+)$/){
					$seq1= $1 . $2;
				}
				#calculate gaps per bp aligned
				my $s1= $seq1;
				my $s2= $seq1;
				#remove all gaps from one string...
				$s1=~ s/-//g;
				my $len1= length($s1);
				#remove all non-gaps from other string...
				$s2=~ s/([A-Za-z])//g;
				my $len2= length($s2);
				my $checklen= length($seq1);
				my $add= $len1 + $len2;
				unless ($checklen == $add){
					print "Error!  Length: $add doesn't match expected $checklen for $spp $locus $sppkey\n";
					die;
				}
				my $gapperbp= ($len2/$len1);
				$gapperbp= sprintf("%.2f", $gapperbp);
				#...and store non-gap length...
				my $spplen= $len1;

				my @avgall;
				my @avgnonself;
				my %sppvals;
				#now, loop through all species
				for my $key2 (sort {$a cmp $b} (keys %seqnames)){
					#print "Comparing to: $key2...\n";
					#create arrays to hold avg. % identity values for comparisons to individual species...
					my @vals1;
					for my $key3 (keys %{$seqnames{$key2}}){
						#if it's the same as 'query' seq...skip
						if ($key3 eq $sppkey){
							next;
						}
						else{
							#print "$key3\n";
							my $seqobj2= $aln-> get_seq_by_id($key3);
							#create a new temp alignment with 2 seqs of interest
							my $tempaln= Bio::SimpleAlign-> new();
							$tempaln-> add_seq($seqobj1);
							$tempaln-> add_seq($seqobj2);
							#remove any gap-only columns...
							$tempaln= $tempaln-> remove_columns(['all_gaps_columns']);
							#get avg. percent identity & store in arrays...
							my $ident= $tempaln-> average_percentage_identity();
							push(@vals1, $ident);
							push(@avgall, $ident);
							if ($key2 ne $sppkey){
								push(@avgnonself, $ident);
							}
						}
					}
					#now that we've looped through all seqs for the current 'target' species...
					#get final 'species-species' average value & store in %sppvals;...
					#first, check array not empty (happens for comparisions to 'self' species when only 1 input seq)
					my $scalval= scalar(@vals1);
					if ($scalval >= 1){
						my $avgident= mean(\@vals1);
						$avgident= sprintf("%.2f", $avgident);
						$sppvals{$key2}= $avgident;
					}
					else{
						$sppvals{$key2}= 'n/a';
					}
				}
				#after looping through all species...calculate 'total' pairwise averages
				my $allident= mean(\@avgall);
				$allident= sprintf("%.2f", $allident);
				my $nonselfident= mean(\@avgnonself);
				$nonselfident= sprintf("%.2f", $nonselfident);

				#now, print out all values for this 'query' seq...
				print OUT "batch$i\t$locus\t$sppkey\t$spplen\t$gapperbp\t$allident\t$nonselfident";
				for my $k1 (@taxa){
					if (exists $sppvals{$k1}){
						print OUT "\t$sppvals{$k1}";
					}
					else{
						print OUT "\tn/a";
					}
				}
				print OUT "\n";
			}
		}
	}
}					
