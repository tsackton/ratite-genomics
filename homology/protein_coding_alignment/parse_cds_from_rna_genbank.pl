#!/usr/bin/perl
use strict;
use warnings;

use Bio::SeqIO;

#create SeqIO input stream to read GenBank file
my $spp= shift(@ARGV);
chomp($spp);
my $infile= "cds_fastas/$spp.rna.gbk";
my $seqio_in= Bio::SeqIO-> new(-file=> $infile, -format=> 'GenBank');

#open SeqIO output stream to write cds
my $outfile= "$spp\_cds.fa";
my $seqio_out= Bio::SeqIO-> new(-file=> ">$outfile", -format=> 'fasta');

my $count= 0;
my $cdscount= 0;
while (my $seqobj= $seqio_in-> next_seq()){
	$count++;
	if ($count=~ /^(\d+)(0{3})$/){
		print "Processed: $count\n";
	}
	my $accession= $seqobj-> accession();
	my $version= $seqobj-> version();
	my $id= $spp . '_' . $accession . '.' . $version;
	my $desc;
	my $string;

	for my $feat_object ($seqobj->get_SeqFeatures) {
		if ($feat_object->primary_tag eq "CDS") {
			$cdscount++;
     			#get the sequence string for this CDS
			#first, try to get the spliced seq [b/c some CDS like anaPla NM_001289841.1 do have]
			$string= uc($feat_object->spliced_seq->seq);
			#if nothing...grab 'seq' from start-end of CDS annotation
			unless(defined($string)){
				$string= uc($feat_object->seq->seq);
			}

			#get the codon start
			if ($feat_object-> has_tag('codon_start')){
				#for my $pos ($feat_object-> get_tag_values('codon_start')){
					#print "codon start: $pos\n";
				#}
				#instead of looping through (as above)...expect only 1...check
				my @pos= $feat_object-> get_tag_values('codon_start');
				my $scalar_pos= scalar(@pos);
				if ($scalar_pos == 1){
					my $codon= shift(@pos);
					#if codon start is 1...nothing to do...
					if ($codon == 1){
						$string= $string;
					}
					#if 2, pad sequence string with 2 Ns
					elsif ($codon == 2){
						$string= 'NN' . $string;
						#print "$id: codon start $codon\n";
					}
					#if 3, pad with 1 N
					elsif ($codon == 3){
						$string= 'N' . $string;
						#print "$id: codon start $codon\n";
					}
					#if anything else...error?
					else{
						print "Error! $id:  Unexpected value for codon start: $codon\n";	
					}
				}
				else{
					print "Error!  $id: Multiple/no values for codon start (@pos)\n";
				}
			}
			else{
				print "$id: no codon start tag...check!\n";
			}
			#NB- if partial gene...'codon start' can be 1 even if translation not in this frame???
			#so, check current string vs. translation
			if ($feat_object-> has_tag('translation')){
				my $checkseqobj= Bio::Seq-> new(-id=> $id, -seq=> $string);
				my $useframe;
				#translate in 1st, 2nd & 3rd forward reading frames...
				for (my $frame= 0; $frame <= 2; $frame++){
					my $prot= $checkseqobj-> translate(-frame => $frame)-> seq;
					#remove terminal stop
					if ($prot=~ /(.+)(\*)$/){
						$prot= $1;
					}
					#substitute internal stops
					$prot=~ s/\*/X/g;

					#remove initial X (from padded Ns at start)
					if ($prot=~ /^(X)(.+)$/){
						$prot= $2;
					}
					for my $translation ($feat_object-> get_tag_values('translation')){
						#20 standard amino acids plus U for selenocysteine and O for pyrrolysine. 
						#The remaining four letters are reserved for ambiguity: B = D or N, J = I or L, X = unknown, Z = E or Q
						#replace nonstandard U & O with X before checking
						$translation=~ s/U/X/g;
						$translation=~ s/O/X/g;

						if ($prot eq $translation){
							$useframe= $frame;
						}
						#in a few cases, we'll have a 3' partial, so 1 extra AA in translation provided in GenBank record
						#relative to translating the extracted CDS region...
						else{
							#chop off last character of GenBank translation & check again...
							chop($translation);
							if ($prot eq $translation){
								$useframe= $frame;
							}
						}	
					}
				}
				unless (defined($useframe)){
					#print a warning
					print "Check: $id (no forward reading frame matches translation\n";
				}
				#if we have frame 1, add 2 Ns to start of seq...
				if ($useframe == 1){
					$string= 'NN' . $string;
				}
				elsif ($useframe == 2){
					$string= 'N' . $string;
				}
			}
			#after padding start of sequence...check if we need to pad end...
			my $triples= sprintf("%.2f", length($string)/3);
			#if 'even'...nothing to add...
			if ($triples=~ /^(\d+)(\.)(0{2})$/){
				$string= $string;
			}
			elsif ($triples=~ /^(\d+)(\.)(33)$/){
				$string= $string . 'NN';
			}
			elsif ($triples=~ /^(\d+)(\.)(67)$/){
				$string= $string . 'N';
			}
			else{
				print "Error!  $id: Unexpected 'triples' value: $triples\n";
			}
			
			#get the gene annotation, or set to 'n/a'
	      		if ($feat_object->has_tag('gene')) {
				my @genes= $feat_object-> get_tag_values('gene');
				my $scalargenes= scalar(@genes);
				if ($scalargenes == 0){
					$desc= "gene='n/a'";
				}
				elsif ($scalargenes == 1){
					my $name= shift(@genes);
					$desc= 'gene=' . $name;
				}
				else{
					print "Error!  $id has $scalargenes gene names (@genes)\n";
				}
	      		}
			my $new_seqobj= Bio::Seq-> new(-id=> $id, -desc=> $desc, -seq=> $string);
			$seqio_out-> write_seq($new_seqobj);
		}
	}
}
print "Parsed: $count GenBank records from $spp\n";
print "Wrote: $cdscount cds seqs\n";
