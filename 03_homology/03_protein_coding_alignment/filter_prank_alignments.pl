#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::AlignIO;
use Bio::SimpleAlign;
use File::Basename;
use File::Copy;

#pass the list of loci
my $list= shift(@ARGV);
chomp($list);

#store the batch & loci to run...
my %files;
open(LOCI, "<$list") || die "Cannot open $list for reading: $!\n";
while(<LOCI>){
	chomp($_);
	my @spla= split(/\t/, $_);
	#store keyed by locus name, with batch as value
	$files{$spla[1]}= $spla[0];
}
my $scal= scalar(keys %files);
print "Stored list of: $scal prank alignments to filter\n";

#logfile for completed jobs on serial requeue stop/restart
my %donejobs;
my $baselist= basename($list);
my $logfile= 'jarvis_dna_filter_logfile_' . $baselist;
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(.+?)(\sdone)$/){
			$donejobs{$1}= 'exists';
		}
	}
	close(DONE);
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

my $processed= 0;
for my $locus (sort {$a cmp $b} (keys %files)){
	my $batch= $files{$locus};
	
	#skip if already done
	if (exists $donejobs{$locus}){
		next;
	}
	#otherwise, run jarvis prank filtering
	else{
		#directories to output filtered alignments (1 for unaltered jarvis filter output, 2nd for filtered plus removal of triplet gap-only cols
		my $topdir1= "/n/regal/edwards_lab/acloutier/oma_duplicate_genes/PAML_filtered_PRANK_fastas/";
		unless (-e $topdir1){
			mkdir $topdir1 || die "Cannot create $topdir1: $!\n";
		}
		my $outdir1= $topdir1 . $batch . '/';
		unless (-e $outdir1){
			mkdir $outdir1 || die "Cannot create $outdir1: $!\n";
		}

		my $topdir2= "/n/regal/edwards_lab/acloutier/oma_duplicate_genes/PAML_final_PRANK_fastas/";
		unless (-e $topdir2){
			mkdir $topdir2 || die "Cannot create $topdir2: $!\n";
		}
		my $outdir2= $topdir2 . $batch . '/';
		unless (-e $outdir2){
			mkdir $outdir2 || die "Cannot create $outdir2: $!\n";
		}

		#print "Processing: $locus\n";
		my $fasta_in= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/PAML_aligned_PRANK_fastas/' . $batch . '/' . $locus . '.fa.best.fas';
		my $filter= `perl /n/regal/edwards_lab/acloutier/oma_duplicate_genes/filter_alignment_fasta_v1.3B.pl --input $fasta_in`;
		#filtering script will output 2 files: 
		#locus.fa.best-15-0.3.filter & locus.fa.best-15-0.3.stat.xls; move these to outdir1
		my $outfile1= "$locus\.fa.best-15-0.3.filter";
		my $outfile2= "$locus\.fa.best-15-0.3.stat.xls";

		my $outfile1b= "$outdir1$locus\.fa.best-15.0.3.filter";
		my $outfile2b= "$outdir1$locus\.fa.best-15.0.3.stat.xls";
		
		#print "Moving jarvis dna filter outfiles to $outdir1\n";
		move($outfile1,$outfile1b) || die "Move failed for $locus filter file: $!\n";
		move($outfile2,$outfile2b) || die "Move failed for $locus stats file: $!\n";

		#now, we'll read in the jarvis filter file & remove wholly uninformative columns,
		#but only in triplets & only in-frame (reason for doing this is to avoid many warning messages in RAxML run for building PAML guide tree)
		my $alnio_in= Bio::AlignIO-> new(-file=> $outfile1b, -format=> 'fasta');
		my %toremove;
		#print "Removing non-informative triplet columns...\n";
		while(my $aln= $alnio_in-> next_aln()){
			#grab alignment slices for each 3 columns
			my $aln_len= $aln-> length();
			for (my $i= 1; $i< $aln_len; $i+=3){
				my $i2= $i + 1;
				my $j= $i + 2;
				#slice current 'window' & use boolean to retain gap-only seqs (otherwise, will throw error if slice contains all gaps for seqs)
				my $slice= $aln-> slice($i,$j,1);
				#print "Checking slice: $i-$j\n";
				
				#get all characters in this alignment slice (except gaps are skipped)
				my @symbolchars= $slice->symbol_chars;
				my $flag= 'no';
				for my $char (@symbolchars){
					$char= uc($char);
					#print "Character: $char\n";
					if (($char ne '-') && ($char ne 'N')){
						$flag= 'yes';
					}
				}
				#only if we didn't have any nongap/N characters...record these columns as ones to remove
				if ($flag eq 'no'){
					$toremove{$i}= 'exists';
					$toremove{$i2}= 'exists';
					$toremove{$j}= 'exists';
				}
			}
			my $newaln= $aln;
			#now, remove columns beginning at END of alignment (& nb 'remove_columns' is 0-based)
			for my $pos (sort {$b <=> $a} (keys %toremove)){
				$pos--;
				$newaln= $newaln-> remove_columns([$pos,$pos]);
			}
			#...and output as SeqIO stream rather than AlignIO to avoid LocatableSeq coords being appended to end of ids
			my $fastaout= $outdir2 . $locus . '.fa';
			my $seqio_out= Bio::SeqIO-> new(-file=> ">$fastaout", -format=> 'fasta');
			for my $seqobj ($newaln-> each_seq()){
				#& make sure all uppercase
				my $idb= $seqobj-> display_id();
				my $seqb= $seqobj-> seq();
				$seqb= uc($seqb);
				#if we have no ACGT letters after filtering...omit sequence
				unless($seqb=~ /([ACGT])/){
					print "Omitting: $idb from $locus (no ACGT bases)\n";
					next;
				}
				my $seqobjb= Bio::Seq-> new(-id=> $idb, -seq=> $seqb, -alphabet=> 'dna');
				$seqio_out-> write_seq($seqobjb);
			}
		}
		#print "$locus done\n";
		print LOG "$locus done\n";
		$processed++;
		if ($processed=~ /^(\d+)(0{2})$/){
			print "Processed: $processed\n";
		}
	}
}
close(LOG);
print "Finished processing: $processed fastas\n";
