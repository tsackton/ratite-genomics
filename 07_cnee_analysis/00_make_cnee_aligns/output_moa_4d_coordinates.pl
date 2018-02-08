#!/usr/bin/perl
use strict;
use warnings;

open(OUT, ">anoDid.4dparsed.bed") || die "Cannot open outfile for writing: $!\n";

my %coords;
#store all of the moa-droNov mapped coordinates, keyed by numerical portion of id
for (my $i= 1; $i <= 23; $i++){
	my $stored= 0;
	my $skipped= 0;
	my $coordfile= 'moa-droNov_coordinates_to_retrieve_' . $i;
	print "Storing info from: $coordfile\n";
	open(COORDS, "<$coordfile") || die "Cannot open $coordfile for reading: $!\n";
	while(<COORDS>){
		chomp($_);
		my @split= split(/\t/, $_);
		my ($id, $scaff, $start, $end, $strand) = ($split[0], $split[1], $split[2], $split[3], $split[4]);
		#numerical portion of id
		my $numid;
		if ($id=~ /^(fourD\.)(\d+)$/){
			$numid= $2;
		}
		else{
			print "Error parsing number from: $id\n";
			die;
		}
		#check that emu base was mapped to a single moa base;
		#otherwise, we won't store & it will be considered missing
		unless ($start == $end){
			$skipped++;
			next;
		}
		else{
			#make current 1-based start be 0-based for .bed format
			$start--;
			my $bedline= "$scaff\t$start\t$end\t$id\t0\t$strand";
			$coords{$numid}= $bedline;
			$stored++;
		}
	}
	close(COORDS);
	print "Stored: $stored coords for file$i (skipped $skipped where emu doesn't correspond to a single moa base)\n";
}
my $totstored= scalar(keys %coords);
print "Stored: $totstored bases total\n";

#for every 4d site...
my $output= 0;
my $havemoa= 0;
my $missing= 0;
for (my $site= 1; $site <= 5739749; $site++){
	#if we have a stored coordinate...print
	if (exists $coords{$site}){
		print OUT "$coords{$site}\n";
		$havemoa++;
		$output++;
	}
	else{
		print OUT "missing\t.\t.\tfourD.$site\t0\t.\n";
		$missing++;
		$output++;
	}
}
print "Output: $output coords\n";
print "Have moa mapped for: $havemoa, missing $missing\n";
