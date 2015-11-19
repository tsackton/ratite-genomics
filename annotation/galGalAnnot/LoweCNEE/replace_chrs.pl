#!/bin/env perl

use strict;
use warnings;

open my $acc, "GCF_000002315.3.assembly.txt" or die;
my %key;
while (<$acc>) {
	chomp;
	next if /^#/;
	my @fields = split(/\t/, $_);
	my $chr = $fields[9];
	$key{$chr} = $fields[6];
	print "$fields[9]\t\t$fields[6]\n"
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
