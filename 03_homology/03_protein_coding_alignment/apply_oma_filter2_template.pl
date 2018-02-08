#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;

#pass the batch name/number (e.g. batch1)...
my $batch= shift(@ARGV);
chomp($batch);

#glob all of the aligned fastas in that batch
my @glob= </n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned_filter1/$batch/*.fa>;

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart)
my $logfile= 'apply_jarvis_filter2_logfile_' . $batch;
my %done;
if (-e $logfile){
	open(DONEJOBS, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONEJOBS>){
		if ($_=~ /^(.+?)(\s)(done)$/){
			$done{$1}= 'exists';
		}
	}
	close(DONEJOBS);

	#reopen for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

my $outdir= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned_filter2/' . $batch . '/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}

my $processed= 0;
for my $filepath (sort {$a cmp $b} (@glob)){
	my $basename= basename($filepath);
	my $locus;
	if ($basename=~ /^(.+)(_aln_filter1)(\.fa)$/){
		$locus= $1;
	}
	elsif ($basename=~ /^(.+)(_filter1)(\.fa)$/){
		$locus= $1;
	}
	else{
		print "Error parsing locus from: $basename\n";
		die;
	}
	$processed++;
	if ($processed=~ /^(\d+)(0{3})$/){
		print "Processed: $processed\n";
	}

	#if we've already processed this locus...skip
	if (exists $done{$locus}){
		next;
	}

	#create an output SeqIO stream to output seqs.
	my $fasta_out= $outdir . $locus . '_aln_filter2.fa';
	my $seqio_out= Bio::SeqIO-> new(-file=> ">$fasta_out", -format=> 'fasta');

	#path to the jarvis 'WrongA.txt' outfile holding info. about any regions to filter...
	my $jpath= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/jarvis_filter_scripts/' . $batch . '/W12S1G6_' . $locus . '_aln_filter1_WrongA.txt';
	#if file doesn't exist...print 'warning' statement (but don't kill script b/c there were some alignments that failed initially & had to be re-run...
	unless (-e $jpath){
		print "Skipping: $locus (no jarvis filter 'WrongA.txt' logfile)\n";
		print LOG "$locus done\n";
		next;
	}
	#otherwise...open for reading & store info. about regions to filter
	my %regions;
	open(JARFILTER, "<$jpath") || die "Cannot open $jpath for reading: $!\n";
	while(<JARFILTER>){
		chomp($_);
		my @split1= split(/\t/, $_);
		my $id1= $split1[1];
		my $start1= $split1[2];
		my $end1= $split1[3];
		#print "Storing: $id1 $start1-$end1\n";
		$regions{$id1}{$start1}= $end1;
	}
	close(JARFILTER);

	#read in the 'filter1' alignment as an AlignIO input stream
	my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
	#create a new alignment object...
	my $alntemp= Bio::SimpleAlign-> new();
	while (my $aln= $alignio_in-> next_aln()){
		if ($aln-> is_flush()){
			#go through seqs. 1-by-1 & mask over any regions to filter...
			foreach my $seqobj ($aln-> each_seq()){
				my $id= $seqobj-> display_id();
				my $seq= $seqobj-> seq();
				
				#if this sequence has bases to filter...
				if (exists $regions{$id}){
					#print "\nSeq: $id\n";
					#explode the filter1 seq on empty string...
					my @bases= split(//, $seq);
					my $count1= 0;
					my %truepos;
					#map each of the positions to its 'true' coordinate in sequence...
					for (my $i= 0; $i <= $#bases; $i++){
						my $pos= $bases[$i];
						if ($pos ne '-'){
							$count1++;
							$truepos{$count1}= $i;
						}
					}
					#now, mask each of the regions to filter...
					for my $startpos (sort {$a <=> $b} (keys %{$regions{$id}})){
						my $endpos= $regions{$id}{$startpos};
						#print "Overwriting: $startpos-$endpos with gaps...\n";
						#get the corresponding array-based coordinates
						my $maskstart= $truepos{$startpos};
						my $maskend= $truepos{$endpos};
						#print "Corresponding array-based coords are: $maskstart-$maskend\n";
						for (my $z= $maskstart; $z <= $maskend; $z++){
							#overwrite the array base with a gap character...
							$bases[$z]= '-';
							#print "Overwriting: $z\n";
						}
					}
					#then, add modified sequence to temp alignment
					my $string= join('', @bases);
					my $new_seqobj= Bio::LocatableSeq-> new(-id=> $id, -seq=> $string);
					$alntemp-> add_seq($new_seqobj);
				}
				#if there were no bases to filter, just add to temp alignment
				else{
					$alntemp-> add_seq($seqobj);
				}
			}
		}
		else{
			print "Error!  Alignment: $locus is not flush!\n";
			die;
		}
	}
	#now, remove any gap-only columns created during filtering...
	$alntemp= $alntemp-> remove_columns(['all_gaps_columns']);
	#...and, output!
	#do as a SeqIO stream rather than AlignIO to avoid LocatableSeq coords appended to end of Seq IDs
	foreach my $seqobj2 ($alntemp-> each_seq()){
		#as long as there are any bases...output; otherwise, omit
		my $checkseq= $seqobj2-> seq();
		my $checkid= $seqobj2-> display_id();
		if ($checkseq=~ /^(-+)$/){
			print "Omitting: $checkid from $locus (no bases remaining)\n";
		}
		else{
			$seqio_out-> write_seq($seqobj2);
		}
	}
	print LOG "$locus done\n";
}
close(LOG);	
print "Finished processing: $processed alignments\n";		 
