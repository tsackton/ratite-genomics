#!/bin/env perl

use strict;
use warnings;

open my $acc, "galGal4.chr2acc" or die;
my %key;
while (<$acc>) {
	chomp;
	next if /^#/;
	my ($chr, $acc) = split;
	$chr = "chr$chr";
	$key{$chr} = $acc;
}

my $bed = shift;
open my $infile, $bed or die;
my $fixed = "$bed.fixed";
open my $outfile, ">$fixed" or die;

while (<$infile>) {
	chomp;
	my @fields = split;
	my $chr = $fields[0];
	my $acc = exists($key{$chr}) ? $key{$chr} : "NA";
	next if $acc eq "NA";
	$fields[0] = $acc;
	print $outfile join("\t", @fields);
	print $outfile "\n";
} 
