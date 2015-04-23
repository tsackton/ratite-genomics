Contiguity
---------

We use the assemblathon_stats.pl (forked here, original from [Keith Bradnam] (https://github.com/ucdavis-bioinformatics/assemblathon2-analysis/blob/master/assemblathon_stats.pl)) script to get genome contiguity measures from FASTA files, as so:

Code to process all fasta or gzipped fasta files in a directory is:

	for FILE in $(ls *.{fa,fas.gz})
	do 
		SP=${FILE%%.*}
    	assemblathon_stats.pl -csv -graph -n 10 $FILE > $SP.out
	done	

Completeness
------------

We use the BUSCO approach for estimating completeness (forked here, original from [Simão et al](http://buscos.ezlab.org/)).

Code to process a single genome is:

	SP=<species name>
	GENOME=/path/to/genome
	BUSCO_v1.0.py -o $SP -in $GENOME -l V -m genome -c 8 -sp chicken
