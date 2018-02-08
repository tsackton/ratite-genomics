#!/usr/bin/perl
use strict;
use warnings;

#pass input NCBI galGal4 annotations in GFF format
my $infile= shift(@ARGV);

#outfile to hold CDS
my $outfile1= 'galGal4_CDS.bed';
open(OUT1, ">$outfile1") || die "Cannot open $outfile1 for writing: $!\n";

#outfile to hold ALL annotated exons
my $outfile2= 'galGal4_all_exons.bed';
open(OUT2, ">$outfile2") || die "Cannot open $outfile2 for writing: $!\n";

#outfile to hold CODING exons
my $outfile3= 'galGal4_coding_exons.bed';
open(OUT3, ">$outfile3") || die "Cannot open $outfile3 for writing: $!\n";

#outfile to hold all introns annotated between coding exons
my $outfile4= 'galGal4_introns.bed';
open(OUT4, ">$outfile4") || die "Cannot open $outfile4 for writing: $!\n";

open(INFILE, "<" . $infile) || die "Cannot open $infile for reading: $!\n";
my %genes;
my %genes2;
my %rnas;
my %mrnas;
my %rnacheck;
my %cds;
my %types;
while(<INFILE>){
	chomp($_);
	if ($_=~ /^(\#)/){
		next;
	}
	my @split= split(/\t/, $_);
	#if this is a 'gene' annotation line...store info in %genes
	if ($split[2] eq 'gene'){
		my $geneid;
		#my $dbxref;
		my $id2;
		my $genename;
		my $biotype;
		#parse the ID string
		if ($_=~ /(ID=)(.+?)(\;)(.+?)(GeneID\:)(.+?)(\;)(Name=)(.+?)(\;)(.+?)(gene_biotype=)(.+?)(\;)/){
			$geneid= $2;
			$id2= $6;
			$genename= $9;
			$biotype= $13;
			$types{$biotype}= 'exists';
			#store the info in %genes, keyed by gene ID (e.g. 'gene1')
			$genes{$geneid}{'id'}= $id2;
			$genes{$geneid}{'name'}= $genename;
			$genes{$geneid}{'biotype'}= $biotype;
		}
		elsif ($_=~ /(ID=)(.+?)(\;)(.+?)(GeneID\:)(.+?)(\;)(Name=)(.+?)(\;)(.+?)(gene_biotype=)(.+?)$/){
			$geneid= $2;
			$id2= $6;
			$genename= $9;
			$biotype= $13;
			$types{$biotype}= 'exists';
			$genes{$geneid}{'id'}= $id2;
			$genes{$geneid}{'name'}= $genename;
			$genes{$geneid}{'biotype'}= $biotype;
		}
		else{
			print "Error parsing: $_\n";
		}
	}
	#if this is an mRNA, store info linking mRNA to parent gene in %genes
	elsif ($split[2] eq 'mRNA'){
		if ($_=~ /(ID=)(.+?)(\;)(Parent=)(.+?)(\;)(.+?)(transcript_id=)([A-Z]{2})(_)(\d+)(\.)(\d)$/){
			my $parent= $5;
			my $rnaid= $2;
			my $transcript_id= $9 . $10 . $11 . $12 . $13;
			
			if (exists $genes{$parent}){
				my $gettype= $genes{$parent}{'biotype'};
				$rnacheck{$rnaid}= $gettype;
				#add a field to store transcript_id for this mRNA linked to parent gene id
				$genes2{$rnaid}{'transcript'}= $transcript_id;
				$genes2{$rnaid}{'parent_gene'}= $parent;

				#if this is an RNA for a protein-coding gene only...
				if ($genes{$parent}{'biotype'} eq 'protein_coding'){
					$rnas{$rnaid}{'parent'}= $parent;
					$rnas{$rnaid}{'id'}= $transcript_id;
				}
			}
			else{
				print "Error!  Parsing a mRNA with no stored 'Parent' info! (Parent: $parent)\n";
				print "$_\n\n";
				die;
			}
		}
		else{
			print "Error parsing: $_\n";
			die;
		}
	}
	#if this is a 'transcript' or a noncoding RNA, etc. don't want to store exons/CDS to retrieve for phylogenies, but do want to record names, etc.
	elsif (($split[2] eq 'transcript') or ($split[2] eq 'ncRNA') or ($split[2] eq 'tRNA') or ($split[2] eq 'V_gene_segment') or ($split[2] eq 'C_gene_segment') or ($split[2] eq 'rRNA')){
		my $rid;
		my $gpar;
		my $tid= '';
		my $ids= $split[8];
		if ($ids=~ /(ID=)(.+?)(\;)/){
			$rid= $2;
		}
		if ($ids=~ /(Parent=)(.+?)(\;)/){
			$gpar= $2;
		}
		if ($ids=~ /(transcript_id=)([A-Z]{2})(_)(\d+)(\.)(\d)/){
			$tid= $2 . $3 . $4 . $5 . $6;
		}
		$genes2{$rid}{'transcript'}= $tid;
		$genes2{$rid}{'parent_gene'}= $gpar;
	}

	#if this is a coding region, 	
	elsif ($split[2] eq 'CDS'){
		my @splitb= split(/\t/, $_);
		my $source= $splitb[0];
		my $start= $splitb[3];
		my $end= $splitb[4];
		#NB- bed starts are 0-based, and ends are 1-based, so decrement GFF start by 1, and leave end as is
		$start--;
		#$end++;
		my $strand= $splitb[6];
		my $phase= $splitb[7];
		my $idstring= $splitb[8];
		#parse the idstring to get the parent mRNA and protein id...
		#if this is a partial CDS, or has frameshift/stop, or is mitochondrial... 
		#it's not a 'good' CDS for our purposes...so just store the count of how many CDS occur for this mRNA, but don't store the CDS itself
		#NB- originally did this including 'Notes', but now leave these ones in (some just have subs. btw. RefSeq & genome...)
		if (($idstring=~ /(partial=true)/) or ($idstring=~ /(transl_except=)/) or ($idstring=~ /(transl_table=2)/) or ($idstring=~ /(exception=rearrangement)/)){
			#get the parent id
			if ($idstring=~ /(Parent=)(.+?)(\;)/){
				my $useparent= $2;
				if (exists $mrnas{$useparent}{'CDS_count'}){
					$mrnas{$useparent}{'CDS_count'}++;
				}
				else{
					$mrnas{$useparent}{'CDS_count'}= 1;
				}
			}
			next;
		}
		if ($idstring=~ /(ID=)(.+?)(\;)(Parent=)(.+?)(\;)(.+?)(protein_id=)([A-Z]{2})(_)(\d+)(\.)(\d)/){
			my $cdsid= $2;
			my $parent2= $5;
			my $protid= $9 . $10 . $11 . $12 . $13;
			if (exists $mrnas{$parent2}{'CDS_count'}){
				$mrnas{$parent2}{'CDS_count'}++;
			}
			else{
				$mrnas{$parent2}{'CDS_count'}= 1;
			}
			
			#store the info. in %cds, keyed first by cds ID, then by start & end (to maintain a unique entry for each individual region)
			#BUT- only do if we have info. stored for a parent mRNA (i.e. these CDS are from a protein-coding gene)
			if (exists $rnas{$parent2}){
				$cds{$cdsid}{$start}{$end}{'parent_rna'}= $parent2;
				#& get the parent gene as well
				my $parentgene= $rnas{$parent2}{'parent'};
				$cds{$cdsid}{$start}{$end}{'parent_gene'}= $parentgene;
				$cds{$cdsid}{$start}{$end}{'source'}= $source;
				$cds{$cdsid}{$start}{$end}{'strand'}= $strand;
				$cds{$cdsid}{$start}{$end}{'id'}= $protid;
				$cds{$cdsid}{$start}{$end}{'phase'}= $phase;
			}
			else{
				my $checktype= $rnacheck{$parent2};
				unless ($checktype eq 'protein_coding'){
					print "No info. stored for: $parent2 (biotype: $checktype)\n";
				}
			}
		}
		else{
			print "Error parsing: $idstring\n";
			die;
		}
	}
	#if this is an exonic region...
	elsif ($split[2] eq 'exon'){
		#we just want to print it back out in .bed format
		my $chr= $split[0];
		my $s= $split[3];
		my $e= $split[4];
		#decrement start coordinates for 0-based .bed (& leave end 1-based)
		$s--;
		
		my $strand2= $split[6];
		#create a name to ID this exon
		#parse the parent RNA ID
		my $string= $split[8];
		if ($string=~ /^(ID=)(.+?)(\;)(Parent=)(.+?)(\;)/){
			my $exonid= $2;
			my $rnaparent= $5;
			
			if (exists $genes2{$rnaparent}){
				my $geneparent= $genes2{$rnaparent}{'parent_gene'};
				my $rna_trans= $genes2{$rnaparent}{'transcript'};
				if (defined($geneparent)){
					my $genepname= $genes{$geneparent}{'name'};
					my $genepid= $genes{$geneparent}{'id'};
					my $bedname;
					if ($rna_trans ne ''){
						$bedname= "Exon=$exonid,gene=$genepname,GeneID=$genepid,transcript_id=$rna_trans";
					}
					else{
						$bedname= "Exon=$exonid,gene=$genepname,GeneID=$genepid";
					}
					print OUT2 "$chr\t$s\t$e\t$bedname\t0\t$strand2\n";
				}
				else{
					if ($rna_trans ne ''){
						print OUT2 "$chr\t$s\t$e\tExon=$exonid,transcript_id=$rna_trans\t0\t$strand2\n";
					}
					else{
						print OUT2 "$chr\t$s\t$e\tExon=$exonid\t0\t$strand2\n";
					}
				}
			}
			#for pseudogenes, no stored RNA parent, so check if there's a stored gene parent instead
			elsif (exists $genes{$rnaparent}){
				my $genepname2= $genes{$rnaparent}{'name'};
				my $genepid2= $genes{$rnaparent}{'id'};
				my $bedname2= "Exon=$exonid,gene=$genepname2,GeneID=$genepid2";
				print OUT2 "$chr\t$s\t$e\t$bedname2\t0\t$strand2\n";
			}
			else{
				print "Error!  No stored parent gene for: $rnaparent (exon: $exonid)\n";
				die;
			}
		}
		else{
			print "Error parsing exon ID: $_\n";
			die;
		}
	}
}
close(INFILE);	

#so, now we have all of the 'good' protein-coding CDS stored...
#we want to output: individual CDS 'exons', concatenated 'cDNA' CDSs, AND introns that separate 'good' CDSs...	
#so...go through all stored CDS...
for my $cds_key (sort {$a cmp $b} (keys %cds)){
	my $scalar_starts= scalar(keys %{$cds{$cds_key}});
	#if there's a single region...no need to calculate introns...so just output this CDS to both 'single' & 'cDNA' files
	if ($scalar_starts == 1){
		for my $start_key (keys %{$cds{$cds_key}}){
			#double-check there's a single end coordinate...
			my $scalar_ends= scalar(keys %{$cds{$cds_key}{$start_key}});
			if ($scalar_ends == 1){
				for my $end_key (keys %{$cds{$cds_key}{$start_key}}){
					#first, double-check that we have the expected # of CDS for this rna parent
					my $p_rna= $cds{$cds_key}{$start_key}{$end_key}{'parent_rna'};
					my $check_parentcds= $mrnas{$p_rna}{'CDS_count'};
					unless (($check_parentcds == $scalar_starts) and ($check_parentcds == $scalar_ends)){
						print "\nError for stored CDS: $cds_key\n";
						print "Stored CDS info: #CDS=$scalar_starts,$scalar_ends\n";
						print "Stored mRNA info: #CDS=$check_parentcds\n";
						die;
					}

					my $p_gene= $cds{$cds_key}{$start_key}{$end_key}{'parent_gene'};
					my $chr_val= $cds{$cds_key}{$start_key}{$end_key}{'source'};
					my $strand_val= $cds{$cds_key}{$start_key}{$end_key}{'strand'};
					my $protid_val= $cds{$cds_key}{$start_key}{$end_key}{'id'};
					my $phase_val= $cds{$cds_key}{$start_key}{$end_key}{'phase'};
					
					#and, get the transcript id...
					my $trans_val= $rnas{$p_rna}{'id'};

					#...and the gene name and id...
					my $p_gene_name= $genes{$p_gene}{'name'};
					my $p_gene_id= $genes{$p_gene}{'id'};

					#double-check if phase 0...
					if ($phase_val != 0){
						print "Warning: Single CDS gene has phase $phase_val (GeneID:$p_gene_id,Name:$p_gene_name)\n";
					}

					my $bed_nameval= "CDS=$cds_key,gene=$p_gene_name,GeneID=$p_gene_id,transcript_id=$trans_val,protein_id=$protid_val,phase=$phase_val";
					print OUT3 "$chr_val\t$start_key\t$end_key\t$bed_nameval\t0\t$strand_val\n";

					my $blocksize= ($end_key - $start_key);
					print OUT1 "$chr_val\t$start_key\t$end_key\t$bed_nameval\t0\t$strand_val\t$start_key\t$start_key\t0\t$scalar_starts\t$blocksize\t0\n";
				}
			}
			else{
				print "Error!  $scalar_ends CDS 'ends' for start $start_key in CDS: $cds_key\n";
			}
		}
	}
	#if we have multiple regions...
	else{
		my @starts;
		my @ends;
		my @blockstarts;
		my @blocksizes;
		my @phases;
		my $p_gene2;
		my $chr_val2;
		my $strand_val2;
		my $protid_val2;
		my $trans_val2;
		my $p_gene_name2;
		my $p_gene_id2;
		my $p_rna2;

		for my $start_key2 (sort {$a <=> $b} (keys %{$cds{$cds_key}})){
			push(@starts, $start_key2);
			#double-check there's a single end coordinate...
			my $scalar_ends2= scalar(keys %{$cds{$cds_key}{$start_key2}});
			if ($scalar_ends2 == 1){
				for my $end_key2 (sort {$a <=> $b} (keys %{$cds{$cds_key}{$start_key2}})){
					push(@ends, $end_key2);

					$p_rna2= $cds{$cds_key}{$start_key2}{$end_key2}{'parent_rna'};
					$p_gene2= $cds{$cds_key}{$start_key2}{$end_key2}{'parent_gene'};
					$chr_val2= $cds{$cds_key}{$start_key2}{$end_key2}{'source'};
					$strand_val2= $cds{$cds_key}{$start_key2}{$end_key2}{'strand'};
					$protid_val2= $cds{$cds_key}{$start_key2}{$end_key2}{'id'};
					my $phase_val2= $cds{$cds_key}{$start_key2}{$end_key2}{'phase'};
					push(@phases, $phase_val2);
					
					#get the transcript id...
					$trans_val2= $rnas{$p_rna2}{'id'};

					#...and the gene name and id...
					$p_gene_name2= $genes{$p_gene2}{'name'};
					$p_gene_id2= $genes{$p_gene2}{'id'};

					#create a 'name' field for output...
					my $bed_nameval2= "CDS=$cds_key,gene=$p_gene_name2,GeneID=$p_gene_id2,transcript_id=$trans_val2,protein_id=$protid_val2,phase=$phase_val2";
					print OUT3 "$chr_val2\t$start_key2\t$end_key2\t$bed_nameval2\t0\t$strand_val2\n";
					#calculate the block size, and the block start relative to first CDS start...
					my $blocksize2= ($end_key2 - $start_key2);
					push(@blocksizes, $blocksize2);
					my $firststart= $starts[0];
					my $blockstart2= ($start_key2 - $firststart);
					push(@blockstarts, $blockstart2);
				}
			}
			else{
				print "Error!  $scalar_ends2 CDS 'ends' for start $start_key2 in CDS: $cds_key\n";
			}
		}
		#now that we've gone through all starts/ends...
		#first, double-check that we have the expected # of CDS for this RNA parent
		my $check_parentcds2= $mrnas{$p_rna2}{'CDS_count'};
		my $checkends= scalar(@ends);
		unless (($check_parentcds2 == $scalar_starts) and ($check_parentcds2 == $checkends)){
			print "\nError for stored CDS: $cds_key\n";
			print "Stored CDS info: #CDS=$scalar_starts,$checkends\n";
			print "Stored mRNA info: #CDS=$check_parentcds2\n";
			die;
		}

		#get the very first start & very last end
		my @starts2= @starts;
		my $chr_start= shift(@starts2);
		my @ends2= @ends;
		my $chr_end= pop(@ends2);

		my $joinstarts= join(',', @blockstarts);
		my $joinsizes= join(',', @blocksizes);
		#double-check that the last start plus last size equals the chromosome end position minus chromosome start
		my $lastblockstart= pop(@blockstarts);
		my $lastblocksize= pop(@blocksizes);
		my $lastsum= $lastblockstart + $lastblocksize;
		my $span= $chr_end - $chr_start;
		unless ($lastsum == $span){
			print "Error!  CDS: $cds_key chromosome 'span' $span (= $chr_end - $chr_start) does not equal last block start plus last block size ($lastblockstart + $lastblocksize = $lastsum)\n";
			die;
		}
		my $firstphase;
		if ($strand_val2 eq '+'){
			$firstphase= shift(@phases);
		}
		elsif ($strand_val2 eq '-'){
			$firstphase= pop(@phases);
		}
		else{
			print "Error!  Unexpected strand: $strand_val2 in CDS: $cds_key\n";
			die;
		}
		my $bed_nameval3= "CDS=$cds_key,gene=$p_gene_name2,GeneID=$p_gene_id2,transcript_id=$trans_val2,protein_id=$protid_val2,phase=$firstphase";
		print OUT1 "$chr_val2\t$chr_start\t$chr_end\t$bed_nameval3\t0\t$strand_val2\t$chr_start\t$chr_start\t0\t$scalar_starts\t$joinsizes\t$joinstarts\n";
		
		#then, figure out intron coordinates between each coding exon
		my $stored_end;
		my $tempval= $scalar_starts;
		for (my $i= 0; $i < $scalar_starts; $i++){
			#if it's the first loop, just store this end
			if ($i == 0){
				$stored_end= $ends[$i];
			}
			else{
				$tempval--;
				#otherwise, get the current start
				my $currstart= $starts[$i];
				#set the intron start on chromosome as the stored end (do not add 1 b/c .bed starts are 0-based, but ends 1-based, so works out...)
				my $usestart= $stored_end;
				#set the intron end on chromosome as CDS start (do not subtract 1 b/c .bed end values are 1-based)
				my $useend= $currstart;
				#if we're on the - strand, want intron numbers to 'increment' in other direction
				my $intnum;
				if ($strand_val2 eq '+'){
					$intnum= $i;
				}
				elsif ($strand_val2 eq '-'){
					$intnum= $tempval;
				}
					
				my $bed_nameval4= "CDS=$cds_key,gene=$p_gene_name2,GeneID=$p_gene_id2,transcript_id=$trans_val2,protein_id=$protid_val2,intron=$intnum";
				print OUT4 "$chr_val2\t$usestart\t$useend\t$bed_nameval4\t0\t$strand_val2\n";
				#then, store this end for next loop...
				$stored_end= $ends[$i];
			}
		}
	}		
}
close(INFILE);
close(OUT1);
close(OUT2);
close(OUT3);
close(OUT4);
