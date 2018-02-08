#!/usr/bin/perl
use strict;
use warnings;

#list of all taxa (excluding chicken reference & moa)
my @taxa= qw( allMis allSin anaPla anoCar aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut );

for my $spp (@taxa){
	print "\nParsing 4d sites for: $spp\n";
	#path to .bed output from liftover of galGal 4d sites
	my $bed_in= $spp . '.4d.bed';
	my $processed= 0;	
	my $output= 0;
	my $missing= 0;
	my $multiple= 0;
	my %multiples;
	my %output_multi;
	my $storednum= 0;
	my $storedline= '';
	open(OUT, ">$spp\.4dparsed.bed") || die "Cannot open output file for writing: $!\n";
	open(BED, "<$bed_in") || die "Cannot open $bed_in for reading: $!\n";

	#loop through a 1st time to record all multiple ids
	my $storeda= 0;
	while(<BED>){
		chomp($_);
		my @splita= split(/\t/, $_);
		my $ida= $splita[3];
		my $numa;
		if ($ida=~ /^(fourD\.)(\d+)$/){
			$numa= $2;
		}
		else{
			print "Error parsing number from: $ida\n";
			die;
		}
		if ($numa == $storeda){
			$multiples{$storeda}= 'exists';
		}
		$storeda= $numa;
	}
	close(BED);
	my $scalara= scalar(keys %multiples);
	print "Stored: $scalara duplicate liftover loci for $spp\n";

	#now, loop through again to print out .bed with 1 entry for each galGal 4d site
	open(BED2, "<$bed_in") || die "Cannot open $bed_in for reading: $!\n";
	while(<BED2>){
		chomp($_);
		$processed++;
		my @split= split(/\t/, $_);
		my $id= $split[3];
		my $num;
		if ($id=~ /^(fourD\.)(\d+)$/){
			$num= $2;
		}
		else{
			print "Error parsing number from: $id\n";
			die;
		}
		
		#if this is the first record...
		if ($processed == 1){
			#if it's not the liftover for 4d site '1', add appropriate amount of 'missing' bases
			if ($num != 1){
				for (my $i= 1; $i < $num; $i++){
					print OUT "missing\t.\t.\tfourD.$i\t0\t.\n";
					$missing++;
					$output++;
				}
			}
			#then, store current num & line for next loop
			$storednum= $num;
			$storedline= $_;
		}
		else{
			#if the current number is stored + 1, means we should output the previous stored line...unless that was a duplicate liftover
			my $plusone= $storednum + 1;
			if ($num == $plusone){
				unless (exists $multiples{$storednum}){
					print OUT "$storedline\n";
					$output++;
				}
			}
			#if current number is same as stored...this is a duplicated liftover...
			elsif ($num == $storednum){
				#only output a line if we don't already have one output for this 4d site
				unless(exists $output_multi{$num}){
					print OUT "multiple\t.\t.\tfourD.$num\t0\t.\n";
					$multiple++;
					$output++;
					$output_multi{$num}= 'exists';
				}
			}
			#if current number is more than stored + 1...output stored line (unless a duplicate), & then appropriate number of 'missing'
			elsif ($num > $plusone){
				unless (exists $multiples{$storednum}){
					print OUT "$storedline\n";
					$output++;
				}
				my $num2;
				for (my $j= $plusone; $j < $num; $j++){
					#print "Printing 'missing' record for site $j\n";
					print OUT "missing\t.\t.\tfourD.$j\t0\t.\n";
					$missing++;
					$output++;
					$num2= $j;
				}
			}
			else{
				print "Error!  Unexpected comparison of current 'num' ($num) to stored/plusone/greater than plusone\n";
				die;
			}
			#then, store current values for next loop
			$storednum= $num;
			$storedline= $_;
		}
	}
	close(BED);
	#if last liftover site wasn't the final galGal site (5739749)...add 'missing' to end
	if ($storednum != 5739749){
		unless (exists $multiples{$storednum}){
			print OUT "$storedline\n";
			$output++;
		}
		print "Processed: $processed .bed lines for $spp...adding 'missing' to end...\n";
		$storednum++;
		for (my $z= $storednum; $z <= 5739749; $z++){
			print OUT "missing\t.\t.\tfourD.$z\t0\t.\n";
			$missing++;
			$output++;
		}
	}
	#if last liftover site was the final galGal site...output, or code as multiple if duplicated
	else{
		unless (exists $multiples{$storednum}){
			print OUT "$storedline\n";
			$output++;
		}
		else{
			print OUT "multiple\t.\t.\tfourD.$storednum\t0\t.\n";
			$multiple++;
			$output++;
		}
	}
	print "Output: $output 4d sites for $spp\n";
	print "(incl. $missing missing sites & $multiple multiple liftover sites [will be recoded as -/N])\n";
	close(OUT);
}		
