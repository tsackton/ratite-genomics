#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::Seq;
use Math::Round;

#lists of species 'groups' (outgroups, neognaths, palaeognaths)
my @list1= qw( allMis anoCar chrPic );
my @list2= qw( anaPla aptFor aquChr balReg calAnn chaPel chaVoc colLiv corBra cucCan egrGar falPer ficAlb fulGla galGal geoFor halLeu lepDis melGal melUnd mesUni nipNip picPub pseHum pygAde serCan taeGut );
my @list3= qw( aptHaa aptOwe aptRow casCas cryCin droNov eudEle notPer rheAme rhePen strCam tinGut );

#map each list to keys of hash
my %spp1 = map { $_ => 1 } @list1;
my %spp2 = map { $_ => 1 } @list2;
my %spp3 = map { $_ => 1 } @list3;

#pass the batch number...
my $batch= shift(@ARGV);
chomp($batch);

#glob all of the aligned fastas in that batch
my @glob= </n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned/$batch/*.fa>;

#open logfile for writing/appending to keep track of any jobs already completed
#(for serial requeue stop/restart)
my $logfile= 'filter_oma_align_logfile_' . $batch;
my %done;
if (-e $logfile){
	open(DONEJOBS, "<$logfile") || die "Cannot open $logfile for reading: $!\n";
	while(<DONEJOBS>){
		if ($_=~ /^(.+?)(\s)(done)$/){
			$done{$1}= 'exists';
		}
	}
	close(DONEJOBS);

	#reopen for appending
	open(LOG, ">>$logfile") || die "Cannot open $logfile for appending: $!\n";
}
else{
	open(LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

my $outdir= '/n/regal/edwards_lab/acloutier/oma_duplicate_genes/fastas_aligned_filter1/' . $batch . '/';
unless (-e $outdir){
	mkdir $outdir || die "Cannot create $outdir: $!\n";
}

for my $filepath (sort {$a cmp $b} (@glob)){
	my $basename= basename($filepath);
	my $locus;
	if ($basename=~ /^(.+)(\.fa)$/){
		$locus= $1;
	}
	else{
		print "Error parsing locus from: $basename\n";
		die;
	}

	#if we've already aligned this locus...skip
	if (exists $done{$locus}){
		next;
	}
	
	my $fasta_out= $outdir . $locus . '_filter1.fa';

	#read in the alignment as an AlignIO input stream
	my $alignio_in= Bio::AlignIO-> new(-file=> $filepath, -format=> 'fasta');
	my %cols;
	my %grps;
	while (my $aln= $alignio_in-> next_aln()){
		if ($aln-> is_flush()){
			my $filtered= $aln;
			my $numseqs= $aln-> num_sequences();
			print "Num. seqs. = $numseqs\n";
			my $alnlen= $aln-> length();
			#go through each seq. & record any non-gap positions in alignment...
			foreach my $seqobj ($aln-> each_seq()){
				my $id= $seqobj-> display_id();
				#parse just the tax ID (without protein ID)
				my $taxid;
				if ($id=~ /^([A-Za-z]{6})(_)/){
					$taxid= $1;
				}
				else{
					print "Error parsing taxon ID from: $id\n";
					die;
				}
				my $seq= $seqobj-> seq();
				#explode on empty string to create array of bases...
				my @bases= split(//, $seq);

				for (my $i= 0; $i < $alnlen; $i++){
					my $base= $bases[$i];
					#if it's not a gap character: 
					if ($base ne '-'){
						#initialize/increment count of non-gap bases at this position...
						if (exists $cols{$i}){	
							$cols{$i}++;
						}
						else{
							$cols{$i}= 1;
						}
						#and also record which 'group' each species belongs to
						if (exists $spp1{$taxid}){
							$grps{$i}{'outgroup'}= 'exists';
						}
						elsif (exists $spp2{$taxid}){
							$grps{$i}{'neognaths'}= 'exists';
						}
						elsif (exists $spp3{$taxid}){
							$grps{$i}{'palaeognaths'}= 'exists';
						}
						else{
							print "Error!  No group membership for: $taxid\n";
							die;
						}
					}
				}
			}
			#now, get list of all columns to remove...
			my %remove;
			#set cutoff at 30% (i.e. need at least 30% of seqs with a (non-gap) base for that column to be retained...)
			#OR- a min. of 10 seqs.
			#OR- seqs. belonging to at least 2 of the 3 'groups' (outgroups, neognaths, palaeognaths)
			my $cutoff= round($numseqs*0.3);
			if ($cutoff > 10){
				$cutoff= 10;
			}
			for my $key1 (keys %cols){
				my $val1= $cols{$key1};
				my $numgrps= scalar(keys %{$grps{$key1}});
				if ($val1 <= $cutoff){
					if ($numgrps < 2){
						$remove{$key1}= 'exists';
					}
				}
			}
			#now, go through alignment & remove those columns...
			#NB- do starting from alignment END, and alignment numbering is 0-based
			for my $key2 (sort {$b <=> $a} keys (%remove)){
				$filtered= $filtered-> remove_columns([$key2,$key2]);
			}
			#now, output the filtered alignment
			#NB- do as a SeqIO stream to avoid coords if output as AlignIO
			my $seqio_out= Bio::SeqIO-> new(-file=> ">$fasta_out", -format=> 'fasta');
			foreach my $seqobj2 ($filtered-> each_seq()){
				$seqio_out-> write_seq($seqobj2);
			}
		}
		else{
			print "Error!  $locus: alignment is not flush!!!\n";
			die;
		}
	}
	print LOG "$locus done\n";	
}
