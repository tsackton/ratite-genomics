#!/usr/bin/perl
use strict;
use warnings;

use Bio::DB::Fasta;
use Bio::Seq;
use Math::Round;

#create a fastadb of all of the moa extracted seqs.
my $moa_cnees= 'anoDid_cnees.fa';
my $fastadb= Bio::DB::Fasta-> new($moa_cnees);

#list of ambiguous DNA codes
my @ambig= qw( Y R W S K M D V H B X N );

open(OUT, ">anoDid_cnee_coverage") || die "Cannot open outfile for writing: $!\n";
print OUT "ID\tgalGal_length\tdroNov_length\tanoDid_length\tanoDid_length_noN\tanoDid_length_no-ambig\tanoDid_coverage_noN (% of galGal)\n";

#store all galGal lengths
my $galbed= 'final_cnees_long.bed';
open(GAL, "<$galbed") || die "Cannot open galGal cnee coordinates for reading: $!\n";
my %gal;
while(<GAL>){
	chomp($_);
	my @split= split(/\t/, $_);
	my ($start, $end, $id)= ($split[1], $split[2], $split[3]);
	my $len= $end - $start;
	$gal{$id}= $len;
}
close(GAL);

#store the emu lengths
open(EMU, "<droNov_cnees_parsed_liftover.bed") || die "Cannot open droNov cnee coordinates for reading: $!\n";
my %emu;
while(<EMU>){
	chomp($_);
	my @split2= split(/\t/, $_);
	my ($start2, $end2, $id2)= ($split2[1], $split2[2], $split2[3]);
	if ($id2=~ /^(ID=)(.+)$/){
		$id2= $2;
	}
	my $len2= $end2 - $start2;
	$emu{$id2}= $len2;
}
close(EMU);

#for each cnee...
my $count= 0;
my $skipped1= 0;
my $skipped2= 0;
my $output= 0;

for my $key1 (sort {$a cmp $b} (keys %gal)){
	$count++;
	if ($count=~ /^(\d+)(0{3})$/){
		print "Processed: $count\n";
	}
	my $gal_len= $gal{$key1};
	#if there's no info. for emu...
	unless (exists $emu{$key1}){
		$skipped1++;
		next;
	}	
	else{
		my $emu_len= $emu{$key1};
		#get the moa sequence for this cnee
		my $moaid= 'anoDid_' . $key1;
		my $seqobj= $fastadb-> get_Seq_by_id($moaid);
		#if there is no moa sequence
		unless(defined($seqobj)){
			$skipped2++;
			next;
		}
		else{
			my $seq= $seqobj-> seq();
			$seq= uc($seq);
			#get the total length
			my $moalen1= length($seq);
			#get the length minus any Ns
			my $seq2= $seq;
			$seq2=~ s/N//g;
			my $moalen2= length($seq2);
			#calculate as % of galGal cnee length
			my $perc= round(($moalen2/$gal_len)*100);
			#get the length minus any Ns & any ambiguous characters
			for my $base(@ambig){
				$seq2=~ s/$base//g;
			}
			my $moalen3= length($seq2);
			print OUT "$key1\t$gal_len\t$emu_len\t$moalen1\t$moalen2\t$moalen3\t$perc\n";
			$output++;
		}
	}
}
close(OUT);
print "Finished processing: $count galGal cnees\n";
print "Skipped: $skipped1 (no droNov coordinates)\n";
print "Skipped: $skipped2 (no anoDid seq)\n";
print "Output results for: $output cnees\n";
