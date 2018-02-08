#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;
use Bio::Seq;
use Bio::LocatableSeq;
use Bio::AlignIO;
use Bio::SimpleAlign;

#list of all taxa (including chicken reference and moa; NB- put reference galGal as 1st sequence)
my @taxa= qw( galGal allMis allSin anaPla anoCar anoDid aptFor aptHaa aptOwe aptRow balReg calAnn casCas chaPel chaVoc cheMyd chrPic colLiv corBra croPor cryCin cucCan droNov eudEle falPer ficAlb fulGla gavGan halLeu lepDis melGal melUnd mesUni nipNip notPer picPub pseHum pygAde rheAme rhePen strCam taeGut tinGut );

#seqio output stream in fasta format (comment/uncomment appropriate line for 4d alignment vs. cnee alignment)
#my $seqio_out= Bio::SeqIO-> new(-file=> ">allspecies_4d_concatenated.fasta", -format=> 'fasta');
my $seqio_out= Bio::SeqIO-> new(-file=> ">allspecies_cnee_concatenated.fasta", -format=> 'fasta');

#and, also print alignment out in phylip format to scan as check...
#my $alnio_out= Bio::AlignIO-> new(-file=> ">allspecies_4d_concatenated.phy", -format=> 'phylip');
my $alnio_out= Bio::AlignIO-> new(-file=> ">allspecies_cnee_concatenated.phy", -format=> 'phylip');

my $aln= Bio::SimpleAlign-> new();
for my $spp (@taxa){
	print "\nProcessing: $spp\n";
	#path to fasta for each species (comment/uncomment for 4d vs. cnees)
	#my $fasta= "$spp\_4d_concatenated.fasta";
	my $fasta= "cnee_concatenated_seqs/$spp\_concatenated_cnees.fasta";

	my $seqio_in= Bio::SeqIO-> new(-file=> $fasta, -format=> 'fasta');
	while(my $seqobj= $seqio_in-> next_seq()){
		#write to output fasta
		$seqio_out-> write_seq($seqobj);
		#add to allspecies alignment (requires locatableseq object)
		my $seq= $seqobj-> seq();
		my $newseqobj= Bio::LocatableSeq-> new(-id=> $spp, -seq=> $seq);
		$aln-> add_seq($newseqobj);
	}
}
#now, double-check alignment
if ($aln-> is_flush()){
	my $len= $aln-> length();
	print "Total alignment length: $len\n";
	my $numseqs= $aln-> num_sequences();
	print "# seqs: $numseqs\n";
}
else{
	print "Error!  Alignment is not flush!\n";
}
#output in phylip format
$alnio_out-> write_aln($aln);
