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
#print "Basedir: $basedir\n";

#check if any results already output...and then open outfile for writing/appending
my %donejobs;
my $outfile= 'alignment_summary_stats';
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
			if ($scalarsplit == 50){
				$donejobs{$spl[1]}= 'exists';
			}
		}
	}
	close(DONE);
	open(OUT, ">>$outfile") || die "Cannot open $outfile for appending: $!\n";
}
else{
	open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
	print OUT "Batch\tLocus\tAlignment Length (bp; with gaps)\tNum. seqs\tAvg. input seq length (no gaps)\tAvg. % identity\tGaps per bp aligned (total)\tGaps per bp aligned (avg. across spp)";
	for my $spp2 (@taxa){
		 print OUT "\t$spp2 #seqs";
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
	print "Processing: $scalarglob fastas in batch$i\n";
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
		#skip if already done
		if (exists $donejobs{$locus}){
			#print "Skipping: $locus (already done)\n";
			next;
		}

		#go through the alignment...
		my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
		my %seqnames;
		my @lens;
		my $gaps= 0;
		my $bases= 0;
		my @avggaps;
		while (my $aln= $alignio_in-> next_aln()){
			if ($aln-> is_flush()){
				my $alnlen= $aln-> length();
				my $numseqs= $aln-> num_sequences();
				my $avgident= $aln-> average_percentage_identity();
				$avgident= sprintf("%.2f", $avgident);

				foreach my $seqobj ($aln-> each_seq()){
					my $id= $seqobj-> display_id();
					#parse out just the taxon name...
					if ($id=~ /^([A-Za-z]{6})(.+)$/){
						my $shortid= $1;
						#initialize/increment count of how many seqs there are for this species...
						if (exists $seqnames{$shortid}){
							$seqnames{$shortid}++;
						}
						else{
							$seqnames{$shortid}= 1;
						}
					}
					else{
						print "Error parsing taxon id from: $id\n";
						die;
					}
					#get the: number of internal gaps, & number of bases...
					my $seq= $seqobj-> seq();
					#replace leading/trailing gaps
					if ($seq=~ /^(-+?)([^-])(.+)$/){
						$seq= $2 . $3;
					}
					if ($seq=~ /^(.+)([^-])(-+)$/){
						$seq= $1 . $2;
					}
					
					my $s1= $seq;
					my $s2= $seq;
					#remove all gaps from one string...
					$s1=~ s/-//g;
					my $len1= length($s1);
					#remove all non-gaps from other string...
					$s2=~ s/([A-Za-z])//g;
					my $len2= length($s2);
					my $checklen= length($seq);
					my $add= $len1 + $len2;
					unless ($checklen == $add){
						print "Error!  Length: $add doesn't match expected $checklen for $locus\n";
						die;
					}
					$gaps+= $len2;
					$bases+= $len1;
					push(@lens, $len1);
					my $gapsper= $len2/$len1;
					push(@avggaps, $gapsper);
				}
			
				my $avglen= mean(\@lens);
				$avglen= sprintf("%.2f", $avglen);

				my $totgaps= $gaps/$bases;
				$totgaps= sprintf("%.2f", $totgaps);
			
				my $grandgaps= mean(\@avggaps);
				$grandgaps= sprintf("%.2f", $grandgaps);

				print OUT "batch$i\t$locus\t$alnlen\t$numseqs\t$avglen\t$avgident\t$totgaps\t$grandgaps";
				for my $spp2 (@taxa){
					if (exists $seqnames{$spp2}){
						print OUT "\t$seqnames{$spp2}";
					}
					else{
						print OUT "\tn/a";
					}
				}
				print OUT "\n";
			}
			else{
				print "Error!  $locus: alignment is not flush!\n";
				die;
			}
		}
	}
}
close(OUT);				
