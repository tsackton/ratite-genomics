#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use File::Basename;


#store info. for each emu (droNov) cnee [from parsed output of liftover from galGal reference]
print "Parsing coordinate info. for droNov...\n";

my $bedfile= 'droNov_cnees_parsed_liftover.bed';
my %ref;
my %starts;
my %ends;
my $tot= 0;
my $stored= 0;
open (BED, "<" . $bedfile) || die "Cannot open $bedfile for reading: $!\n";
while(<BED>){
	chomp($_);
	$tot++;
	my @split= split(/\t/, $_);
	
	my ($scaff, $start, $end, $id, $strand) = ($split[0], $split[1], $split[2], $split[3], $split[5]);
	if ($id=~ /^(ID=)(.+?)$/){
		$id= $2;
	}

	#increment .bed coordinate for 1-based start
	$start++;

	#store info. in %ref, keyed first by scaffold, then by start, then end, then locus then with strand as value
	$ref{$scaff}{$start}{$end}{$id}= $strand;
	#also store the start coordinate & end in separate hashes..
	$starts{$scaff}{$start}= 'exists';
	$ends{$scaff}{$end}= 'exists';
	$stored++;
}
close(BED);
print "Finished processing: $tot seqs. for droNov\n";
print "Stored: $stored seqs. to map to moa coords\n\n";

#now, map the stored reference regions to their corresponding sequence in the moa genome
print "\nMapping droNov reference-> moa coordinates...\n";

my $outfile= 'moa-droNov_coordinates_to_retrieve';
open(OUT, ">$outfile") || die "Cannot open $outfile for writing: $!\n";
print OUT "Locus\tMoa_scaffold\tStart(1-based)\tEnd(1-based)\tStrand\n";
my %map;
my $scalarscaffs= scalar(keys %ref);
my $s_count= 0;
for my $key1 (sort {$a cmp $b} (keys %ref)){
	$s_count++;
	print "Processing: $key1 ($s_count of $scalarscaffs)\n";
	#create temp hashes to hold last moa coord mapping to reference start and first coord mapping to end for just this scaffold...
	my %starttemp;
	my %endtemp; 
	#get the file holding reference-> moa coordinate mapping
	my $coordfile= '/n/regal/edwards_lab/acloutier/moa_mapping/moa-droNov_bowtie_remapping/individual_scaffold_coordinates_final/' . $key1;
	#if we didn't have any moa reads mapped to this scaffold...
	unless (-e $coordfile){
		#don't throw an error, b/c could be true that we just don't have any mapped moa bases to this reference scaffold
		print "Skipping: $key1 (no reads mapped)\n";
		next;
	}
	else{
		my $storedmoa= 0;
		open(COORDS, "<" . $coordfile) || die "Cannot open $coordfile for reading: $!\n";
		while(<COORDS>){
			chomp($_);
			if ($_=~ /^(Original)/){
				next;
			}
			else{
				my @spl= split(/\t/, $_);
				my $ref= $spl[0];
				my $moa= $spl[2];	
				#if current reference base is stored as a reference genome start...
				if (exists $starts{$key1}{$ref}){
					#store the correspoding moa base in %starttemp, keyed by ref base
					#NB- will overwrite previous entries if multiple occurrences for a given ref base- this is what we want, to keep the last one
					#(e.g. omit any insertions directly before ref base...)
					$starttemp{$ref}= $moa;
					#BUT- if the stored moa base from previous loop is the same as current...means there was a moa insertion or gap 
					#extending to this site...so increment moa coord to use by 1
					if (($moa == $storedmoa) && ($storedmoa != 0)){
						my $moanew= $moa + 1;
						#overwrite stored value from above...
						$starttemp{$ref}= $moanew;
					}
				}
				#if current base is stored as a reference genome end...
				if (exists $ends{$key1}{$ref}){
					#we want to store the FIRST occurrence in %endtemp (e.g. ignore any insertions after this base)...
					if (exists $endtemp{$ref}){
						if ($moa > $endtemp{$ref}){
							next;
						}
						else{
							$endtemp{$ref}= $moa;
						}
					}
					else{
						$endtemp{$ref}= $moa;
					}
				}
				#and, store current moa value for next loop...
				$storedmoa= $moa;
			}
		}
		close(COORDS);
		#now that we've grabbed all coords for this scaffold...store info in %cmap
		#NB- store instead of just printing so that we can print all out in sorted order at end...
		for my $key2 (keys %{$ref{$key1}}){
			#get the moa coord that corresponds to this reference genome start
			my $moastart= $starttemp{$key2};
			for my $key3 (keys %{$ref{$key1}{$key2}}){
				#get the moa coord that corresponds to this reference genome end
				my $moaend= $endtemp{$key3};
				#'key4' is the locus
				for my $key4 (keys %{$ref{$key1}{$key2}{$key3}}){
					#get the strand
					my $moastrand= $ref{$key1}{$key2}{$key3}{$key4};

					#store moa info. in %map, keyed by locus
					$map{$key4}{'scaffold'}= $key1;
					$map{$key4}{'start'}= $moastart;
					$map{$key4}{'end'}= $moaend;
					$map{$key4}{'strand'}= $moastrand;
				}
			}
		}
	}
}

#now, output...
print "\nWriting final mapping coordinates...\n";	
my $outcount= 0;
#sorted by locus name...
for my $k1 (sort {$a cmp $b} (keys %map)){
	my $printscaff= $map{$k1}{'scaffold'};
	my $printstart= $map{$k1}{'start'};
	#if 'start' is 0, make 1
	if ($printstart == 0){
		$printstart= 1;
	}
	my $printend= $map{$k1}{'end'};
	my $printstrand= $map{$k1}{'strand'};
	print OUT "$k1\t$printscaff\t$printstart\t$printend\t$printstrand\n";
	$outcount++;
}
print "Finished writing coordinates for: $outcount loci\n\n";			
close(OUT);
