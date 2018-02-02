#!/usr/bin/perl
use strict;
use warnings;

use POSIX;
use File::Basename;

#get list of input phylip files from command line
my $list= shift(@ARGV);

my $basename= basename($list);

#store the numerical portion of each file
my %nums;
open(LIST, "<$list") || die "Cannot open $list for reading: $!\n";
while(<LIST>){
	chomp($_);
	if ($_=~ /^(cnees_concatenated_)(\d+)(\.phylip)$/){
		$nums{$2}= 'exists';
	}
	else{
		print "Error parsing phylip 'number' from: $_\n";
		die;
	}
}
close(LIST);

#store all 'done' jobs so won't redo on Odyssey stop/restart
#NB- comment/uncomment lines according to which constraint tree is being used
my $logfile= $basename . 'branchlengths_ver1_logfile';
#my $logfile= $basename . 'branchlengths_ver2_logfile';
#my $logfile= $basename . 'branchlengths_ver3_logfile';
my %donejobs;
if (-e $logfile){
	open(DONE, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONE>){
		chomp($_);
		if ($_=~ /^(.+)(\s)(done)$/){
			$donejobs{$1}= 'exists';
		}
	}
	close(DONE);
	#reopen logfile for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

for my $key1 (sort {$a <=> $b} (keys %nums)){
	#skip if already done
	if (exists $donejobs{$key1}){
		next;
	}
	else{
		#output file name
		my $outname= "cnee_branchlengths_ver1_$key1";
		#my $outname= "cnee_branchlengths_ver2_$key1";
		#my $outname= "cnee_branchlengths_ver3_$key1";

		#path to input partitions file
		my $partfile= 'RAxML_branchlength_input/' . $key1 . '_partitions';

		#path to input alignment file
		my $alnfile= 'RAxML_branchlength_input/cnees_concatenated_' . $key1 . '.phylip';

		my $model= 'GTRGAMMA';

		my $treefile= 'ratiteTree_withMoa.ver1.nh';
		#my $treefile= 'ratiteTree_withMoa.ver2.nh';
		#my $treefile= 'ratiteTree_withMoa.ver3.nh';

		#generate random number
		my $range= 100000;
		my $min= 1;
		my $rand1= int(rand($range)) + $min;

		#output directory to hold results for this batch
		my $outdir= 'cnee_branchlengths_ver1_raw/' . $key1 . '/';
		#my $outdir= 'cnee_branchlengths_ver2_raw/' . $key1 . '/';
		#my $outdir= 'cnee_branchlengths_ver3_raw/' . $key1 . '/';
		#clean up any files from partial runs due to Odyssey stop/restart
		if (-e $outdir){
			my $cmd= "rm -r $outdir";
			system($cmd);
		}
		#then...create/recreate directory...
		mkdir ($outdir) || die "Cannot create directory $outdir: $!\n";

		#command to run RAxML
		my $raxcmd= "raxmlHPC-SSE3 -p $rand1 -m $model -f e -M -s $alnfile -q $partfile -t $treefile -n $outname -w $outdir";
		system($raxcmd);

		print LOG "$key1 done\n";
	}
}	
