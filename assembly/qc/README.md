Contiguity
---------

We use the assemblathon_stats.pl (forked here, original from [Keith Bradnam] (https://github.com/ucdavis-bioinformatics/assemblathon2-analysis/blob/master/assemblathon_stats.pl)) script to get genome contiguity measures from FASTA files, as so:

Code to process all fasta or gzipped fasta files in a directory is:

```bash
	for FILE in $(ls *.{fa,fas.gz})
	do 
		SP=${FILE%%.*}
    	assemblathon_stats.pl -csv -graph -n 10 $FILE > $SP.out
	done	
```

Completeness
------------

We use the BUSCO approach for estimating completeness (forked here, original from [Simão et al](http://buscos.ezlab.org/)).

Code to process a batch of genomes in a list is:

```bash
	for SPFA in $(cat <GENOME LIST>)
	 	SP=${SPFA%%.*}
		mkdir $SP
		cd $SP
		ln -s ../metazoa metazoa
		GENOME=/n/regal/edwards_lab/ratites/qc/genomes/seq/$SP.fa.gz
		gzip -c -d $GENOME > $SP.genome
		BUSCO_v1.0.py -o $SP -in $SP.genome -l M -m genome -c 1 -sp chicken --flank 30000
		cd ..
	done
```
**Note:** The BUSCO analysis is ongoing as we revise the input gene set to account for various issues