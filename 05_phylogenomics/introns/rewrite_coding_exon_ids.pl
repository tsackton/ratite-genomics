#!/usr/bin/perl
use strict;
use warnings;

#PURPOSE: rewrite the IDstring in galGal4 coding exons .bed file to include the chicken start/stop coordinates

my $infile= 'galGal4_coding_exons.bed';
open(INFILE, "<" . $infile) || die "Cannot open $infile for reading: $!\n";

my $outfile= 'galGal4_coding_exons_indexedID.bed';
open (OUT, ">" . $outfile) || die "Cannot open $outfile for writing: $!\n";

while(<INFILE>){
	chomp($_);
	my @split= split(/\t/, $_);
	#add the start & end coordinates into the idstring
	if ($split[3]=~ /^(CDS=cds)(\d+?)(,)(.+)$/){
		my $newstring= $1 . $2 . '_' . $split[1] . '-' . $split[2] . $3 . $4;
		#replace old idstring
		$split[3]= $newstring;
		print OUT join("\t", @split);
		print OUT "\n";
	}
	else{
		print "Error parsing ID string: $split[3]\n";
		die;
	}
}
close(INFILE);
close(OUT);
