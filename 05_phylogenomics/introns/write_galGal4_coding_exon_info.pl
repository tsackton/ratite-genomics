#!/usr/bin/perl
use strict;
use warnings;

my $infile= 'input_annotations/galGal4_coding_exons_indexedID.bed';
open(INFILE, "<" . $infile) || die "Cannot open $infile for reading: $!\n";
my %exons;
my %strands;

while(<INFILE>){
	chomp($_);
	my @split= split(/\t/, $_);
	my ($scaffold, $start, $end, $strand, $id)= ($split[0], $split[1], $split[2], $split[5], $split[3]);
	
	my $cds;
	my $gene;
	my $geneid;
	my $transcript;
	my $protein;
	my $phase;
	if ($id=~ /^(CDS=)(cds)(\d+?)(_)(.+?)(,)(gene=)(.+?)(,)(GeneID=)(.+?)(,)(transcript_id=)(.+?)(,)(protein_id=)(.+?)(,)(phase=)(\d)$/){
		$cds= $2 . $3;
		$gene= $8;
		$geneid= $11;
		$transcript= $14;
		$protein= $17;
		$phase= $20;
	}
	else{
		print "Error parsing: $id\n";
		die;
	}

	$exons{$cds}{$start}{'scaffold'}= $scaffold;
	$exons{$cds}{$start}{'gene'}= $gene;
	$exons{$cds}{$start}{'gene_id'}= $geneid;
	$exons{$cds}{$start}{'transcript'}= $transcript;
	$exons{$cds}{$start}{'protein'}= $protein;
	$exons{$cds}{$start}{'end'}= $end;
	$exons{$cds}{$start}{'strand'}= $strand;
	$exons{$cds}{$start}{'phase'}= $phase;

	$strands{$cds}= $strand;
}
close(INFILE);

#go through all exons for each CDS & determine exon rank (e.g. rank in CDS, but there could be noncoding exons outside CDS)
for my $key1 (sort {$a cmp $b} (keys %exons)){
	my $count= 0;
	#get the strand orientation for this cds
	my $dir= $strands{$key1};
	my @ordered;
	#if on the plus strand, want to count exons from lowest -> highest start...
	if ($dir eq '+'){
		@ordered= sort {$a <=> $b} (keys %{$exons{$key1}});
	}
	#if on the minus strand, count highest -> lowest start...
	elsif ($dir eq '-'){
		@ordered= sort {$b <=> $a} (keys %{$exons{$key1}});
	}
	else{
		print "Error!  Unexpected strand: $dir in $key1\n";
		die;
	}

	for my $key2 (@ordered){
		#increment exon count, then store exon rank in %exons
		$count++;
		$exons{$key1}{$key2}{'rank'}= $count;
	}
}

#now, print parsed info. to outfile
open(OUT, ">galGal4_coding_exons_info") || die "Cannot open outfile for writing: $!\n";
print OUT "CDS\tGene\tGeneID\tTranscript\tProtein\tScaffold\tStart\tEnd\tStrand\tPhase\tRank (in CDS)\n";

for my $k1 (sort {$a cmp $b} (keys %exons)){
	for my $k2 (sort {$a <=> $b} (keys %{$exons{$k1}})){
		print OUT "$k1\t$exons{$k1}{$k2}{'gene'}\t$exons{$k1}{$k2}{'gene_id'}\t$exons{$k1}{$k2}{'transcript'}\t$exons{$k1}{$k2}{'protein'}\t$exons{$k1}{$k2}{'scaffold'}\t$k2\t$exons{$k1}{$k2}{'end'}\t$exons{$k1}{$k2}{'strand'}\t$exons{$k1}{$k2}{'phase'}\t$exons{$k1}{$k2}{'rank'}\n";
	}
}	
close(OUT);
