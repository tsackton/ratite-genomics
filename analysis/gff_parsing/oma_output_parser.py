#!/usr/bin/env python

import sys, re, os, string, exceptions, os.path, argparse
from collections import defaultdict

#better option processing
parser = argparse.ArgumentParser(description="Read a directory of OMA output and a directory of parsed GFF files and extract missing proteins.")
parser.add_argument('hog', help='HOG Directory', action='store')
parser.add_argument('gff', help="GFF table output directory", action='store')
parser.add_argument('-v', '--verbose', help='Display status messages', action='store_true') 
parser.add_argument('-p', '--prefix', help="Prefix to use for output files (prefix.assigned and prefix.unassigned)", action='store', default="oma_parse")
parser.add_argument('-g', '--gff_ext', help="Extension of the table files in the GFF dir, defaults to tab", action="store", default="tab")
parser.add_argument('-f', '--fasta_ext', help="Extension of the fasta files in the HOG dir, defaults to fa", action="store", default="fa")
args=parser.parse_args()

input_dir = args.hog
gff_dir = args.gff

###PROCESS GFF TABLES####

#dict with prot to gene
prot_to_gene={}
#dict with gene to species
gene_to_sp={} 
#dict with longest protein for each gene
#if there is more than one the choice is random(-ish), based on the first one we encounter in the field
gene_to_longest={}

for tab in os.listdir(gff_dir):
	if tab.endswith(args.gff_ext):
		gff_to_open = gff_dir + '/' + tab
		with open (gff_to_open) as f:
			for line in f:
				line = line.rstrip('\n')
		
				#skip comment lines
				if line.startswith("#"):
					continue
			
				fields=line.split("\t")
				#fields[2] is gene id
				#fields[4] is biotype
				#fields[8] is cds_len
				#fields[9] is protein id
				#fields[10] is is_longest
			
				if fields[4] == "protein_coding":
					#this is all we care about
					#need to strip versions from protein id because of HOG defline issues
					try:
						prot_id = re.findall('(\w+\-*\w*)\.*\d*', fields[9])[0]
						#does not check for key existence so protein ids with 2 gene ids will get clobbered
						prot_to_gene[prot_id] = fields[2]
					except:
						sys.stderr.write ("Problem getting protein id for " + line + "  in " + gff_to_open + "\n")
					
					gene_to_sp[fields[2]] = fields[0]		
					
					#fill in gene to longest dict
					if (fields[10] == "Y") and (not fields[2] in gene_to_longest):
						#longest protein and key doesn't already exist
						gene_to_longest[fields[2]]=fields[9]							

###############

assigned = args.prefix + '.assigned'
unassigned = args.prefix + '.unassigned'

#input dir to process

#this is a dict that will store everything that is assigned (geneIDs) and get a count of the number of transcripts
hogassigned = defaultdict(list)

#this dict stores prot -> hog information so that we can keep track of whether isoforms are assigned to different HOGs
prot_to_hog={}

for hog in os.listdir(input_dir):
	if hog.endswith(args.fasta_ext):
		hog_to_open = input_dir + '/' + hog

		#all we care about are parsing deflines, so don't bother with the overhead of biopython

		with open (hog_to_open) as f:
			for line in f:
				if line.startswith(">"):
					#fasta defline, parse
					#there are two possible patterns
					#one is >sp|gene|prot [sp]
					#other is >prot.ver [sp]
					
					try:
						#see if we can get a prot id with the simpler thing
						cur_prot=re.findall('>(\w+)\.\d+\s+.*', line)[0]
					except:
						try:
							cur_prot=re.findall('>\w+[|]+\w*[|]+(\w+\-*\w*)\s+.*', line)[0]
						except:
							sys.stderr.write ("Problem getting accession for " + hog_to_open + " at " + line)
							continue
					
					prot_to_hog[cur_prot]=hog

					try:					
						cur_gene = prot_to_gene[cur_prot]
						hogassigned[cur_gene].append(cur_prot)
						if len(hogassigned[cur_gene]) > 1:
							#if there is more than one protein from a gene
							hogcheck=defaultdict(int)
							for prot_to_check in hogassigned[cur_gene]:
								prot_hog_val = prot_to_hog[prot_to_check]
								hogcheck[prot_hog_val] += 1
							if len(hogcheck) > 1:
								sys.stderr.write (cur_gene + " has multiple isoforms in different HOGS: " + "\t".join(hogcheck.keys()) + "\n")
						
					except KeyError:
						sys.stderr.write ("Problem getting gene for " + cur_prot + " at " + hog_to_open + " line: " + line)
					
with open(assigned, 'w') as asg, open(unassigned, 'w') as un:
	for gene in gene_to_sp:
		if gene in hogassigned:
			asg.write(gene + "\t" + gene_to_sp[gene] + "\t" + ",".join(hogassigned[gene]) + "\n")
		else:
			un.write(gene + "\t" + gene_to_sp[gene] + "\t" + gene_to_longest[gene] + "\n")

