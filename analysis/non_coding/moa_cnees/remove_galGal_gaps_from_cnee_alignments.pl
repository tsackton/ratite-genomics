#!/usr/bin/perl
use strict;
use warnings;

use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;

#get the batch number from command line
my $batch= shift(@ARGV);

#path to input directory of aligned cnees
my $indir= "/n/regal/edwards_lab/acloutier/moa_cnees/input_fastas_allspecies_cnees_aligned/$batch/";

#output directory...
my $topdir= 'input_fastas_allspecies_cnees_aligned_no-galGal-gaps/';
unless (-e $topdir){
	mkdir ($topdir) || die "Cannot create $topdir: $!\n";
}
my $outdir= $topdir . $batch . '/';
unless (-e $outdir){
	mkdir($outdir) || die "Cannot create $outdir: $!\n";
}

#store any already completed jobs
my $logfile= '/n/regal/edwards_lab/acloutier/moa_cnees/remove_galGal_gaps_log_' . $batch;
my %donejobs;
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(.+?)(\s)(done)$/){
			$donejobs{$1}= 'exists';
		}
	}
	close(DONE);
	#reopen file for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}
#delete any output fastas that aren't stored in %donejobs
print "\nCleaning up files from Odyssey stop/restart...\n";
my @glob= </n/regal/edwards_lab/acloutier/moa_cnees/$outdir/*.fasta>;
for my $filepath (sort {$a cmp $b} (@glob)){
	my $basename= basename($filepath);
	my $rootname;
	if ($basename=~ /^(.+)(\.)(fasta)$/){
		$rootname= $1;
	}
	else{
		print "Error parsing 'root' name from: $basename\n";
		die;
	}
	unless (exists $donejobs{$rootname}){
		unlink($filepath);
	}
}

#glob all input fastas in batch for processing...
my @glob2= <$indir*.fasta>;
my $scalarglob2= scalar(@glob2);
print "Preparing to process: $scalarglob2 fastas in $batch\n";

my $processed= 0;
for my $filepath2 (sort {$a cmp $b} (@glob2)){
	if ($processed=~ /^(\d+)(0{2})$/){
		print "Processed: $processed\n";
	}
	my $basename2= basename($filepath2);
	my $locus;
	if ($basename2=~ /^(.*)(_aln\.fasta)$/){
		$locus= $1;
	}
	else{
		print "Error parsing locus from: $basename2\n";
		die;
	}
	if (exists $donejobs{$locus}){
		$processed++;
		next;
	}
	else{
		my $alnio_in= Bio::AlignIO-> new(-file=> $filepath2, -format=> 'fasta');
		while (my $aln= $alnio_in-> next_aln()){
			if ($aln-> is_flush()){
				my $aln_len= $aln-> length();

				#hash to hold alignment columns where chicken has a gap character
				#(NB-store as 0-based values, b/c this is what 'remove_columns' method used later will require)
				my %galgaps;

				#get the galGal seq
				my $seqobj= $aln-> get_seq_by_id('galGal');
				my $seq= $seqobj-> seq();
				#split on empty string to explode into array of bases
				my @bases= split(//, $seq);
				for (my $i= 0; $i< $aln_len; $i++){
					my $base= $bases[$i];
					if ($base eq '-'){
						$galgaps{$i}= 'exists';
					}
				}
				my $toremove= scalar(keys %galgaps);
				#now, remove each of these alignment columns, starting from END
				for my $key1 (sort {$b <=> $a} (keys %galgaps)){
					$aln= $aln-> remove_columns([$key1,$key1]);
				}
				my $newlen= $aln-> length();
				my $diff= ($aln_len - $newlen);
				unless ($diff == $toremove){
					print "Error!  $locus: should remove $toremove columns with galGal gap\n";
					print "Orig. aln length= $aln_len, New aln length= $newlen (diff= $diff, not $toremove)\n";
					die;
				}
				#now, output as a SeqIO stream (to avoid BioLocatableSeq coordinates being appended to end of ids)
				#& also make all uppercase...
				my $outfile= $outdir . $locus . '.fasta';
				my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');
				foreach my $seqobj2 ($aln-> each_seq()){
					my $id2= $seqobj2-> display_id();
					my $seq2= $seqobj2-> seq();
					$seq2= uc($seq2);
					my $seqobj3= Bio::Seq-> new(-id=> $id2, -seq=> $seq2, -alphabet=> 'dna');
					$seqio_out-> write_seq($seqobj3);
				}
				print LOG "$locus done\n";
				$processed++;
			}
			else{
				print "Error!  $locus input alignment is not flush!\n";
				die;
			}
		}
	}
}
print "Finished processing: $processed alignments in $batch\n";
				
