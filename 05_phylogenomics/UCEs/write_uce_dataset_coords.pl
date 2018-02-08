#!/usr/bin/perl
use strict;
use warnings;

use Statistics::Basic qw(:all);

my @taxa= qw( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );
my $scalarspp= scalar(@taxa);

#store the full set of UCEs from the galGal4 chicken reference...
my $galbed= 'galGal_uces.bed';
my %galref;
open (GALBED, "<$galbed") || die "Cannot open $galbed for reading: $!\n";
while(<GALBED>){
	chomp($_);
	my @spl= split(/\t/, $_);
	my $glen= $spl[2] - $spl[1];
	if ($spl[3]=~ /^(galGal,)(.+)$/){
		$galref{$2}= $glen;
	}
	else{
		print "Error parsing galGal id: $spl[3]\n";
		die;
	}
}
close(GALBED);
my $gscal= scalar(keys %galref);
print "\nFinished storing lengths for: $gscal galGal reference UCEs\n\n";

#store the 'good' liftovers for each taxon (unique target region across all reference liftovers)
my %loci;
for my $spp (@taxa){
	print "\nProcessing: $spp\n";
	my $bedfile= "halLiftover/$spp\_uces_parsed_liftover.bed";
	unless (-e $bedfile){
		print "Error!  No file: $bedfile\n";
		die;
	}
	open(BEDIN, "<$bedfile") || die "Cannot open $bedfile for reading: $!\n";
	my $count= 0;
	while(<BEDIN>){
		chomp($_);
		my @split= split(/\t/, $_);
		my $len= $split[2] - $split[1];
		if ($split[3]=~ /^(ID=)(.+)$/){
			$loci{$2}{$spp}= $len;
			$count++;
		}
		else{
			print "Error parsing id: $split[3]\n";
			die;
		}
	}
	close(BEDIN);
	print "Stored info for: $count loci\n";
}

my $goodloc= 0;
my %tokeep;
#allow up to 1 missing species/locus...
$scalarspp= $scalarspp - 1;
for my $key1 (sort {$a cmp $b} (keys %loci)){
	my $meanlen;
	my $scalartax= scalar(keys %{$loci{$key1}});
	#&& also check if we have a good ostrich seq (b/c will be tested as outgroup in some analyses...)
	if (($scalartax >= $scalarspp) && (exists $loci{$key1}{'strCam'})){
		$goodloc++;
		$tokeep{$key1}= 'exists';
	}
}
print "Good loci (0 or 1 missing spp.): $goodloc\n";
my $scalar2= scalar(keys %tokeep);
print "Check: $scalar2\n";

#print a list of all of the uce loci to use
open(OUTLIST, ">uces_3158_loci_list") || die "Cannot open outfile for writing: $!\n";
for my $cle (sort {$a cmp $b} (keys %tokeep)){
	print OUTLIST "$cle\n";
}
close(OUTLIST);

#now, for each species output .bed lines for loci we want to keep...
for my $spp2 (@taxa){
	my $bedfile2= "halLiftover/$spp2\_uces_parsed_liftover.bed";
	my $outfile= "halLiftover/$spp2\_uces_to_use.bed";
	open(BEDIN2, "<$bedfile2") || die "Cannot open $bedfile2 for reading: $!\n";
	open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
	
	while(<BEDIN2>){
		chomp($_);
		my @split2= split(/\t/, $_);
		if ($split2[3]=~ /^(ID=)(.+)$/){
			my $id2= $2;
			if (exists $tokeep{$id2}){
				print OUT "$_\n";
			}
		}
		else{
			print "Error parsing: $split2[3]\n";
			die;
		}
	}
	close(BEDIN2);
	close(OUT);
}
#and, do chicken...
open (GALBED2, "<$galbed") || die "Cannot open $galbed for reading: $!\n";
my $outfile2= "halLiftover/galGal_uces_to_use.bed";
open(OUT2, ">$outfile2") || die "Cannot open $outfile2 for writing: $1\n";

while(<GALBED2>){
	chomp($_);
	my @spl2= split(/\t/, $_);
	if ($spl2[3]=~ /^(galGal,)(.+)$/){
		my $loc= $2;
		if (exists $tokeep{$loc}){
			my $newid= 'ID=' . $loc;
			print OUT2 "$spl2[0]\t$spl2[1]\t$spl2[2]\t$newid\t$spl2[4]\t$spl2[5]\n";
		}
	}
}
close(GALBED2);
close(OUT2);
