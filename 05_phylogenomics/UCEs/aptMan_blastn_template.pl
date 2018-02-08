#!/usr/bin/perl
use strict;
use warnings;

#list Bioperl modules used
use Bio::SeqIO;
use Bio::Seq;
use Bio::Tools::Run::StandAloneBlastPlus;
use Bio::Search::Result::BlastResult;
use Bio::SearchIO;
use Bio::Search::SearchUtils;
use Bio::Search::Hit::GenericHit;
use File::Basename;

#get first reference, then query species from command line
my $refspp= shift(@ARGV);

my $queryspp= shift(@ARGV);

#create a StandAloneBlastPlus object
my $database= 'blast_dbs/' . $refspp;

#print "Creating BLAST factory...\n";	
my $factory= Bio::Tools::Run::StandAloneBlastPlus -> new(-db_name=> $database);

#get the input fasta file from the command line
my $fastafile= shift(@ARGV);

#for stops/restarts on serial requeue partition
my $logfile= $refspp . '-' . $queryspp . '_logfile';
#if we've already run some...
my %donejobs;
if (-e $logfile){
	#store all 'done' loci...
	open(DONEJOBS, "<" . $logfile) || die "Cannot open $logfile for reading: $!\n";
	while(<DONEJOBS>){
		chomp($_);
		if ($_=~ /^(.+?)(\s)(done)$/){
			$donejobs{$1}= 'exists';
		}
	}
	close(DONEJOBS);
	#then, reopen for APPENDING
	open (LOG, ">>$logfile") || die "Cannot open $logfile for writing: $!\n";
}
else{
	open (LOG, ">$logfile") || die "Cannot open $logfile for writing: $!\n";
}

#create file to hold text output of blast results
my $file= "$refspp-$queryspp\_blastn_results";
if (-e $file){
	open (RESULTS, ">>" . $file) || die "Cannot open output file to hold blastn results: $!\n";
}
else{
	open (RESULTS, ">" . $file) || die "Cannot open output file to hold blastn results: $!\n";
	print RESULTS "Query\tLength (bp)\tNum hits\tHit name\tAccession\tHit length (bp)\tHit Strand\tEvalue\tNum HSPs\tAligned query length (bp)\tFraction aligned query\tFraction aligned hit\tFraction identical\tGaps\tQuery Start\tQuery End\tHit Start\tHit End\n";
}

#create a SeqIO stream to read through query seqs
my $seqio_in= Bio::SeqIO-> new(-file=> $fastafile, -format=> 'fasta');

my $query_count= 0;
while (my $seqobj= $seqio_in-> next_seq()){
	#parse the locus name from id
	my $locname;
	my $qid= $seqobj-> display_id();
	$query_count++;
	if (exists $donejobs{$locname}){
		next;
	}

	#set a flag to record how many hits we've output (allow a max. of 10)
	my $hitflag= 0;
	my $query_length= $seqobj->length();
			
	#set blast parameters (blastn, min. evalue cutoff) to return a reference to a blast object containing a blast report
	#NB-use NCBI blastn 'somewhat similar' search parameters
	my $result= $factory-> blastn(-query=> $seqobj, -num_alignments=> 100, -method_args=> ['-evalue'=> 1e-10, '-ungapped'=> 0, 'perc_identity'=> 10, '-dbsize'=> 0, 'searchsp'=> 0, '-dust'=> 'yes', '-window_size'=> 0, '-penalty'=> -3, '-reward'=> 2, '-gapopen'=> 5, '-gapextend'=> 2, '-strand'=> 'both', '-word_size'=> 11]);

	my $hitcount= $result->num_hits;
		
	#if there are no hits for current query, skip the remainder of processing
	if (($hitcount eq '-') or ($hitcount == 0)){
		print RESULTS "$locname\t$query_length\tNO HITS\n";
		print LOG "$locname done\n";
		$factory-> cleanup;
		next;
	}
	#otherwise, iterate through each hit & process info.
	else{
		while (my $hit= $result-> next_hit()){
			$hitflag++;
			#just print out the top hit...
			if ($hitflag == 1){
				my $hit_name= $hit-> name();
				my $hit_accession= $hit-> accession();
				my $hit_length= $hit-> length();
				my $num_hsps= $hit-> num_hsps();
				my $hit_strand= $hit-> strand('hit');
				my $evalue= $hit-> significance();
				my $query_start= $hit-> start('query');
				my $query_end= $hit-> end('query');
				my $hit_start= $hit-> start('hit');
				my $hit_end= $hit-> end('hit');
				my $query_aln= $hit-> length_aln('query');
				my $frac_aln= $hit-> frac_aligned_query();
				my $frac_aln_hit= $hit-> frac_aligned_hit();
				my $frac_identical= $hit-> frac_identical('query');
				my $gaps= $hit-> gaps('total');

				print RESULTS "$locname\t$query_length\t$hitcount\t$hit_name\t$hit_accession\t$hit_length\t$hit_strand\t$evalue\t$num_hsps\t$query_aln\t$frac_aln\t$frac_aln_hit\t$frac_identical\t$gaps\t$query_start\t$query_end\t$hit_start\t$hit_end\n";			
			}
			else{
				$factory-> cleanup;
			}
		}
		$factory-> cleanup;
	}
	print LOG "$locname done\n";
}
close(RESULTS);
close(LOG);
print "Finished blasting: $query_count query seqs\n";
