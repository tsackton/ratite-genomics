#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;


#open an output file to hold best tree line from each mpest file
my $output_file= 'mpest_parsed_output_run1_nexus.tre';
open(OUT, ">$output_file") or die "Cannot open $output_file to write results parsed from mpest output files: $!\n";

#open a 2nd output file to write trees in newick format
my $output_file2= 'mpest_parsed_output_run1_newick.tre';
open(OUT2, ">$output_file2") or die "Cannot open $output_file2 to write results parsed from mpest output files: $!\n";

#use perl glob to capture all mpest output files
my @mpest_outfiles= <mpest_output_run1/*.tre>;
my $scalar_mpest= scalar(@mpest_outfiles);
print "\nPreparing to parse $scalar_mpest mpest files\n\n";

my $processed_files= 0;
#iterate through each mpest input file
for my $mpest_file (sort {$a cmp $b} @mpest_outfiles){
	$processed_files++;
	
	my $basename= basename($mpest_file);
	my $filenum;
	if ($basename =~ /^(genetree)(\d+)(\.tre\.tre)$/){
		$filenum= $2;
	}
	else{
		print "Error parsing file # from: $basename\n";
	}

	#if this is the very first file processed, print out the NEXUS header lines
	if ($processed_files == 1){
		print OUT "#NEXUS\n";
		print OUT "Begin trees;\n";
		print OUT "  translate\n";
	}
	
	my %tree_likes;
	#open current file by attaching a filehandle
	open(MPEST, "<" . $mpest_file) or die "Cannot open $mpest_file for reading: $!\n";

	#set a 'flag' variable to check whether we've successfully grabbed the desired line from each file...
	my $flag= 'no';
	#read through file line-by-line
	my %taxa;
	my $taxcount= 0;
	my $index= 0;
	while(<MPEST>){
		chomp($_);
		#if current line holds the 'mpest tree'...
		if ($_ =~ /^(\s\s)(tree)(\s)(mpest)/){
			$index++;
			
			#parse out the likelihood score for this tree...
			my $like;
			if ($_=~ /^(\s\s)(tree)(\s)(mpest)(\s)(\[)(-)(\d+)(\.)(\d+)(\])/){
				$like= $7 . $8 . $9 . $10;
				#substitute the file number for 'mpest' (i.e. to make a valid nexus output file)
				$_ =~ s/mpest/$filenum/;
				#then, store entire tree line, keyed by likelihood score, then by index in %tree_likes
				$tree_likes{$like}{$index}{$_}= 'exists';
			}
			else{
				print "genetree: $filenum Error parsing likelihood score from tree string:\n";
				print "$_\n";
				die;
			}
		}
		else{
			#if we're on the first file processed & current line holds a code-> taxon name...
			if (($processed_files == 1) and ($_=~ /^(\t)(\d+)(\s)([A-za-z]+)(,|;)$/)){
				$taxcount++;
				$taxa{$taxcount}= $_;
			}
		}
	}
	close(MPEST);
	#print out the taxa if first file processed...
	if ($processed_files == 1){
		for my $taxkey (sort {$a <=> $b} (keys %taxa)){
			print OUT "$taxa{$taxkey}\n";
		}
	}
		
	#get the highest scoring tree to print to outfile
	my @likes= sort {$a <=> $b} (keys %tree_likes);
	my $usekey= pop(@likes);	

	#check if there's only 1 tree with this value...if so, print
	my $scalar_best= scalar(keys %{$tree_likes{$usekey}});
	if ($scalar_best == 1){
		for my $bestindex (keys %{$tree_likes{$usekey}}){
			for my $besttree (keys %{$tree_likes{$usekey}{$bestindex}}){
				print OUT "$besttree\n";
				#parse this treeline to output in newick format...
				if ($besttree=~ /^(\s\s)(tree)(\s)(\d+)(\s)(\[)(-)(\d+)(\.)(\d+)(\])(\s)(=)(\s)(.+)/){
					print OUT2 "$15\n";
				}
				else{
					print "Error parsing newick sting from: $besttree\n";
					die;
				}
			}
		}
	}
	#if we had multiple 'best' trees...use the one from last line in mpest genetree file, but print warning that there are multiple trees w same max. like
	else{
		my @sorted= sort {$a <=> $b} (keys %{$tree_likes{$usekey}});
		my $bestindex2= pop(@sorted);
		for my $besttree2 (keys %{$tree_likes{$usekey}{$bestindex2}}){
			print OUT "$besttree2\n";
			if ($besttree2=~ /^(\s\s)(tree)(\s)(\d+)(\s)(\[)(-)(\d+)(\.)(\d+)(\])(\s)(=)(\s)(.+)/){
				print OUT2 "$15\n";
			}
			else{
				print "Error parsing newick sting from: $besttree2\n";
				die;
			}
		}
		print "Warning: multiple 'best' trees with likelihood $usekey found for genetree $filenum...check!!!\n";
	}
}
print OUT "end\;\n";
close(OUT);
print "Finished parsing 'best' mpest trees from $processed_files mpest files...\n";
print "Parsed results written to outfile: $output_file\n\n";
