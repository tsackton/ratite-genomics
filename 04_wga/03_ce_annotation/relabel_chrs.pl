#!/usr/bin/perl

my $key = shift;
my $bed = shift;

open KEY, $key or die;
my $key;
while (<KEY>) {
	chomp;
	my ($ncbi, undef, undef, $ucsc) = split;
	$key{$ncbi} = $ucsc;
}

open BED, $bed or die;
while (<BED>) {
	chomp;
	my @fields = split;
	$fields[0] = $key{$fields[0]};
	my $out = join("\t", @fields);
	print "$out\n";
}
