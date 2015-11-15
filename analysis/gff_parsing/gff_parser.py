#!/usr/bin/env python

import sys, re, os, string, exceptions, itertools, pprint, argparse
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter
import pickle

class Transcript(object):
	def __init__(self, geneID, proteinID):
		self.geneID = geneID
		self.proteinID = proteinID

	def __hash__(self):
		return hash((self.geneID, self.proteinID))

	def __eq__(self, other):
		return (self.geneID, self.proteinID) == (other.geneID, other.proteinID)

class Gene(dict):
	#gene 
	def __init__(self, id):
		self.id=id
		self.acc=None
		self.transcripts=()
		self.biotype=None
		self.name=None
		
	def add_trans(self, transcript):
		self.transcripts.append(transcript)
		
	
	
#better option processing
parser = argparse.ArgumentParser(description="Read a GFF file and parses protein coding genes.")
parser.add_argument('gff', help='GFF file to process', action='store')
parser.add_argument('output', help='Name of output file (which will contain transcript list with gene and length information', action='store')
parser.add_argument('-e', '--exceptions', help='File to store exception data to (default=None)', action='store', default='/dev/null')
parser.add_argument('-p', '--pickle', help='Name of pickle object to store parsed dictionary (default=None)', action='store', default='/dev/null')
parser.add_argument('-s', '--species', help='Species name to use (default=GFF file base name)', action='store')
parser.add_argument('-b', '--biotype', help='Extract biotype information for each gene and store in this file', action='store', default='/dev/null')
args=parser.parse_args()

input = args.gff
output = args.output
error = args.exceptions
picklename = args.pickle
species = args.species
biotype_file = args.biotype
if species == None:
	species = re.findall('(.*?).gff', input)[0]

#done processing options

dictionary = OrderedDict()
#genes is a dictionary of gene classes, keyed on the ID field from the GFF
genes = {}
trans_ids = []

sys.stderr.write ("Reading in " + input + " GFF file\n")

with open (input, 'r') as f, open (error, 'w') as e:
	for line in f:
		#this will currently only work on NCBI GFFs that have ID=rna(\d+) and ID=cds(\d+) in them
		#need to make it flexible to work on any GFF using the codes in the 3rd column
		
		#skip comment lines
		if line.startswith("#"):
			continue
		
		#get gfftype and feature id which should be present for every line
		try:
			gfftype=line.split("\t")[2]
		except:
			e.write ('* NO_GFF_TYPE' + line)
			continue
		
		try:
			#this doesn't work, gets too much in the \S+
			#replace with a split on ;?
			feature_id=re.findall('ID=(\S+);', line)[0]
		except:
			e.write ('* NO_FEATURE_ID' + line)	
			continue
				
		if gfftype == "gene":
			genes[feature_id]=Gene(feature_id)
		
			#get biotype
			biotype=re.findall('gene_biotype=(\w+)[;]*.*$', line)
			if biotype:
				genes[feature_id].biotype=biotype[0]
			else:
				genes[feature_id].biotype="unknown"
				
			try:
				acc=re.findall('(GeneID:\d+)[,|;]', line)[0]
			except:
				acc=feature_id
			
			genes[feature_id].acc=acc
			
			try:
				genename=re.findall('Name=(\w+)[;]*.*$', line)[0]
			except:
				genename=acc
			
			genes[feature_id].name=genename
	
		if gfftype == "mRNA":
			#this should be all protein coding transcripts, except those on mt in some cases (need to deal with this case later)
			#now we need to get an id for this transcript
			if not re.findall('(transcript_id=\S+_\d+.?\d+)', line):
				e.write ('* TRANS ' + line)
				continue
			else:
				RNAID = re.findall('ID=rna(\d+)', line)[0]
				transcript_id = re.findall('(transcript_id=\S+_\d+.?\d+)', line)[0]
				trans = (RNAID, transcript_id)
				trans_ids.append(trans)

		if 'ID=cds' in line:
			if not re.findall('(GeneID:\d+)[,|;]', line):
				e.write ('* GENE ' + line)
				continue
			else:
				geneID = re.findall('(GeneID:\d+)[,|;]', line)[0]

			if not re.findall('(protein_id=\S+_\d+.?\d+)', line):
				e.write ('* PROT ' + line)
				continue
			else:
				proteinID = re.findall('(protein_id=\S+_\d+.?\d+)', line)[0]

			if not re.findall('Parent=rna(\d+)', line):
				e.write ('* PARENTID ' + line)
				continue
			else:
				parent = re.findall('Parent=rna(\d+)', line)[0]

			nums = re.findall('CDS\s(\d+)\s(\d+)', line)
			start = int(nums[0][0])
			end = int(nums[0][1])
			length = end - start

			if Transcript(geneID, proteinID) in dictionary:
				dictionary[Transcript(geneID, proteinID)][0] += length
			else:
				dictionary[Transcript(geneID, proteinID)] = [length, 'none']
				dictionary[Transcript(geneID, proteinID)][1] = parent

sys.stderr.write ("Done reading GFF file.\n")
sys.stderr.write ("Writing biotype file, if requested.\n")

#make gene key file
with open(biotype_file, 'w') as b:
	for gffid in genes.keys():
		b.write (genes[gffid].id + "\t" + genes[gffid].acc + "\t" + genes[gffid].name + "\t" + genes[gffid].biotype + "\n")

sys.stderr.write ("Done biotype file, if requested.\n")

splice = []

for thing in dictionary:
	for ids in trans_ids:
		if dictionary[thing][1] == ids[0]:
			dictionary[thing][1] = ids[1]
			break

for key, multi in groupby(dictionary.iteritems(), lambda x: x[0].geneID):
	splice.append(list(multi))

gpDict = OrderedDict()

with open(output, 'w') as f:
	for thing in splice:
		longest = max(thing, key=itemgetter(1))
		for each in thing:
			searchableProteinID = re.findall('protein_id=(\S+_\d+).?\d+', each[0].proteinID)[0]
			gpDict[searchableProteinID] = each[0].geneID
			f.write (species)
			f.write ('\t' + each[0].geneID)
			f.write ('\t' + each[0].proteinID)
			f.write ('\t' + each[1][1])
			f.write ('\t' + str(each[1][0]))
			if len(thing) > 1 and each is longest:
				f.write ('\t*`')
			if len(thing) == 1:
				f.write ('\t*')
			f.write ('\n')
		f.write ('\n')

pickle.dump(gpDict, open(picklename, "wb"))
