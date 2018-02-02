#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::DB::Fasta;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;

my @taxa= qw( galGal allMis allSin anaPla anoCar anoDid aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut );

#get current list from command line 
#(cnee_alignment_batches split into lists of ~126 loci with bash split utility)
my $batch= shift(@ARGV);

my $fastafile= "cnees_concatenated_$batch\.fasta";
my $seqio_out= Bio::SeqIO-> new(-file=> ">$fastafile", -format=> 'fasta');

#store loci for current list
my %batches;
my $listfile= 'concat_lists/list' . $batch;
open(BATCHES, "<$listfile") || die"Cannot open list of batches for reading:$!\n";
while(<BATCHES>){
	chomp($_);
	if ($_=~ /^(Batch)/){
		next;
	}
	else{
		my @spl= split(/\t/, $_);
		#store keyed by batch, then locus
		$batches{$spl[0]}{$spl[1]}= 'exists';
	}
}
close(BATCHES);

#create alignment object to hold concatenated seqs. for each species
my $aln= Bio::SimpleAlign-> new();

for my $spp (@taxa){
	print "\nProcessing: $spp\n";
	#store all missing & duplicated liftovers for this species
	my %missing;
	my %multiple;
	my $parseinfo;
	if ($spp eq 'anoDid'){
		$parseinfo= 'moa_final_cnees_long_liftover_parsing_log';
	}
	else{
		$parseinfo= 'final_cnees_long_liftover_parsing_log';
	}
	open(PARSE, "<$parseinfo") || die "Cannot open $parseinfo for reading: $!\n";
	while(<PARSE>){
		chomp($_);
		if ($_=~ /^(Locus)/){
			next;
		}
		else{
			my @spl2= split(/\t/, $_);
			my ($id2, $spp2, $comment2) = ($spl2[0], $spl2[1], $spl2[2]);
			if ($spp2 eq $spp){
				if ($comment2 eq 'target_region_short'){
					next;
				}
				elsif ($comment2 eq 'no_liftover'){
					$missing{$id2}= 'exists';
				}
				elsif (($comment2 eq 'target_region_too_long') || ($comment2 eq 'multiple_liftover_regions')){
					$multiple{$id2}= 'exists';
				}
				else{
					print "Error!  Unexpected parsing comment: $comment2 for $spp2 $id2\n";
					die;
				}
			}
		}
	}
	close(PARSE);
	my $num_missing= scalar(keys %missing);
	my $num_multiple= scalar(keys %multiple);
	#print "Stored: $num_missing missing loci & $num_multiple duplicated/overly long loci for $spp\n";

	#string to hold concatenated output for this species
	my $string= '';

	#go through each locus in each batch...in same order for all species so we can concatenate seqs across species after...
	my $retrieved= 0;
	my $addmissing= 0;
	my $addmultiple= 0;
	my $count= 0;

	for my $batchkey (sort {$a cmp $b} (keys %batches)){
		for my $lockey (sort {$a cmp $b} (keys %{$batches{$batchkey}})){
			$count++;
			#path to alignment for this locus
			my $alnfile= "input_fastas_allspecies_cnees_aligned_no-galGal-gaps/$batchkey/$lockey\.fasta";
			my $alnio_in= Bio::AlignIO-> new(-file=> $alnfile, -format=> 'fasta');
			while (my $aln= $alnio_in-> next_aln()){
				if ($aln-> is_flush()){
					my $seq;
					my $aln_len= $aln-> length();
					#try to get the sequence for target species
					my $seqobj= $aln-> get_seq_by_id($spp);
					if (defined($seqobj)){
						$seq= $seqobj-> seq();
						$retrieved++;
					}
					else{
						#if this is a missing seq...append a string of gap characters of equal length to alignment
						if (exists $missing{$lockey}){
							$seq= '-' x $aln_len;
							$addmissing++;
							#print "Missing: $lockey\n";
						}
						#if multiple liftover regions, or region too long...use string of Ns
						elsif (exists $multiple{$lockey}){
							$seq= 'N' x $aln_len;
							$addmultiple++;
							#print "Multiple: $lockey\n";
						}
						else{
							print "Error!  No alignment seq for $lockey, but not stored as missing or multiple!!!\n";
							#NB-could be true for moa (anoDid) b/c very short scaffolds don't appear in final assembly
							if ($spp eq 'anoDid'){
								#so, put as missing
								#print "Writing $lockey as 'missing'\n";
								$seq= '-' x $aln_len;
								$addmissing++;
							}
						}
					}
					$seq= uc($seq);
					$string= $string . $seq;
				}
				else{
					print "Error!  Alignment $lockey is not flush!\n";
					die;
				}
			}
		}
	}
	my $seqobj_out= Bio::Seq-> new(-id=> $spp, -seq=> $string, -alphabet=> 'dna');
	$seqio_out-> write_seq($seqobj_out);
	my $seqlen= length($string);
	print "$spp length: $seqlen\n";
}
#now, read in the fasta file & output in phylip format
my $alnio_in= Bio::AlignIO-> new(-file=> $fastafile, -format=> 'fasta');
my $phylip= "cnees_concatenated_$batch\.phylip";
my $alnio_out= Bio::AlignIO-> new(-file=> ">$phylip", -format=> 'phylip');
while (my $alnobj= $alnio_in-> next_aln()){
	if ($alnobj-> is_flush()){
		$alnio_out-> write_aln($alnobj);
	}
}
#& delete unnecessary fasta...
unlink($fastafile);
