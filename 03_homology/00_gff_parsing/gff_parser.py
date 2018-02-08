#!/usr/bin/env python

import sys, re, os, string, exceptions, itertools, pprint, argparse, pickle, gzip

class Transcript(dict):
	#transcript
	def __init__(self, id, parent, type):
		self.parent = parent
		self.type = type
		self.acc = None
		self.id = id
		self.len = 0
		self.cdslen = 0
		self.protid = None
		
		
class Gene(dict):
	#gene 
	def __init__(self, id):
		self.id=id
		self.acc=None
		self.transcripts={}
		self.biotype=None
		self.name=None
		self.protlengths=()
		self.tranlengths=()
			
	def add_trans(self, trans_id, transcript, overwrite=None):
		if not trans_id in self.transcripts:
			self.transcripts[trans_id]=transcript
		elif overwrite:
			self.transcripts[trans_id]=transcript
		
	def get_trans(self, trans_id):
		return self.transcripts[trans_id]
		
	def update_trans(self, trans_id, attr, val):
		#no error checking!
		setattr(self.transcripts[trans_id], attr, val)
			
	def add_length(self, trans_id, len, type):
		if type == "CDS":
			self.transcripts[trans_id].cdslen += len
		elif type == "exon":
			self.transcripts[trans_id].len += len
			
		if self.transcripts[trans_id].len < self.transcripts[trans_id].cdslen:
			self.transcripts[trans_id].len = self.transcripts[trans_id].cdslen
	
	def is_longest_prot(self, trans_id):
		if self.transcripts[trans_id].cdslen == 0:
			return "<NA>"
			
		if not self.protlengths:
			self.protlengths = [self.transcripts[x].cdslen for x in self.transcripts.keys()]
		
		if self.transcripts[trans_id].cdslen >= max(self.protlengths):
			return "Y"
		else:
			return "N"

def xstr(s):
	if s is None:
		return "<NA>"
	else:
		return str(s)

#better option processing
parser = argparse.ArgumentParser(description="Read a GFF file and parses protein coding genes.")
parser.add_argument('gff', help='GFF file to process', action='store')
parser.add_argument('-t', '--table', help='File to store table of information from GFF', default='/dev/null', action='store')
parser.add_argument('-e', '--exceptions', help='File to store exception data to (default=None)', action='store', default='/dev/null')
parser.add_argument('-p', '--pickle', help='Name of pickle object to store parsed dictionary (default=None)', action='store', default='/dev/null')
parser.add_argument('-s', '--species', help='Species name to use (default=GFF file base name)', action='store')
parser.add_argument('-v', '--verbose', help='Display status messages', action='store_true') 
args=parser.parse_args()

#get file names
input_file = args.gff
table_file = args.table
error_file = args.exceptions
pickle_file = args.pickle

#get other variables
species = args.species

if species == None:
	species = re.findall('(.+)\.gff', input_file)[0]

#done processing options

#genes is a dictionary of gene classes, keyed on the ID field from the GFF
genes = {}
#trans to parent links transcript feature ids to their gene id parent for ease of lookup
trans_to_parent = {}

#transcript types holds the gff features that are considered transcripts
#note that this includes C_gene_segment and V_gene_segment, as they are equivalent to transcripts except for rearrangements

transcript_types = frozenset(["mRNA", "ncRNA", "rRNA", "tRNA", "transcript", "C_gene_segment", "D_gene_segment", "J_gene_segment", "V_gene_segment", "primary_transcript"])

if args.verbose: 
	sys.stderr.write ("Reading in " + input_file + " GFF file\n")

with open (input_file, 'r') as f, open (error_file, 'w') as e:
	for line in f:
		
		line = line.rstrip('\n')
		
		#skip comment lines
		if line.startswith("#"):
			continue
		
		#get list of gff fields, verify 8 fields / tab-delimited
		try:
			gff_fields=line.split("\t")
		except:
			e.write ('* PARSE_FAILURE ' + line + "\n")
			continue
		
		if len(gff_fields) != 9:
			e.write ('* WRONG_LEN ' + line + "\n")
			continue
			
		#get feature id from the 8th field of the list (gff_fields[8])
		#WARNING: does not handle quoting properly
		try:
			gff_features = dict(x.split('=') for x in gff_fields[8].strip(';').split(';'))
		except:
			e.write ('* FEATURE_PARSE_ERR ' + line + "\n")
			continue
		
		try:
			#get ID
			feature_id=gff_features['ID']
		except:
			e.write ('* NO_FEATURE_ID ' + line + "\n")	
			continue
				
		if gff_fields[2] == "gene":
			#create a gene entry in the genes dict and add information
			
			if not feature_id in genes:		
				genes[feature_id]=Gene(feature_id)
		
			try:
				if "pseudo" in gff_features:
					genes[feature_id].biotype = "pseudogene"
				else:
					genes[feature_id].biotype = gff_features['gene_biotype']			
			except:
					genes[feature_id].biotype = None
				
			try:
				genes[feature_id].acc = re.findall('GeneID:(\d+)', gff_features['Dbxref'])[0]
			except:
				genes[feature_id].acc=feature_id
						
			try:
				genes[feature_id].name=gff_features['Name']
			except:
				genes[feature_id].name=genes[feature_id].acc
				
		if gff_fields[2] in transcript_types:
			#we've found a transcript, now we need to get the parent and the id (should already by in feature_id)
			
			try:
				parent_id = gff_features['Parent']
			except:
				if gff_fields[2] in ("tRNA", "rRNA"):
					parent_id = feature_id
				else:
					e.write('* NO_PARENT_ID ' + line + "\n")
					continue
			
			if parent_id in genes:
				#gene exists in genes hash
				genes[parent_id].add_trans(feature_id, Transcript(feature_id, parent_id, gff_fields[2]))
			else:
				#gene does not exist, have to add it
				#this means that the mRNA line comes before the gene line which is odd
				genes[parent_id]=Gene(parent_id)
				genes[parent_id].add_trans(feature_id, Transcript(feature_id, parent_id, gff_fields[2]))
			
			trans_to_parent[feature_id] = parent_id
				
			#now we fill out features of transcript
			#first we will try to get the accession from either the transcript_id field or the name field
			#falling back on feature ID if neither of those is populated
			try:
				genes[parent_id].update_trans(feature_id, "acc", gff_features['transcript_id'])
			except:
				try:
					genes[parent_id].update_trans(feature_id, "acc", gff_features['Name'])
				except:
					genes[parent_id].update_trans(feature_id, "acc", feature_id)
			
			#now set gene biotype of parent if we can figure it out (i.e., if it is none and the gff_fields[2] is not mRNA
			if (genes[parent_id].biotype is None) or (genes[parent_id].biotype == "other"):
				try:
					genes[parent_id].biotype = gff_features['gbkey']
				except:
					if gff_fields[2] == "mRNA":
						genes[parent_id].biotype = "protein_coding"
					elif gff_fields[2] == "ncRNA":
						genes[parent_id].biotype = "misc_RNA"
					elif gff_fields[2] == "tRNA":
						genes[parent_id].biotype = "tRNA"
					elif gff_fields[2] == "rRNA":
						genes[parent_id].biotype = "rRNA"

		if gff_fields[2] in ("CDS", "exon"):
			#now we've arrived at the base level feature
			#first thing to do is get the transcript id for the parent
			#but in some rare cases the parent may not be a transcript, it may be straight to a gene
			
			try:
				parent_id_base = gff_features['Parent']
			except:
				e.write('* NO_PARENT_ID ' + line + "\n")
				continue
			
			elem_len = int(gff_fields[4]) - int(gff_fields[3]) + 1	
			
			try:
				elem_id = gff_features['protein_id']
			except:
				try:
					elem_id = gff_features['transcript_id']
				except:
					try:
						elem_id = gff_features['Name']
					except:
						try: 
							elem_id = gff_features['product']
						except:
							if feature_id.endswith(':cds'):
								elem_id = feature_id[:-4]
							else:
								elem_id = feature_id
			
			#there can be more than one parent for a given exon / CDS, so need to convert parent_id to a list
			#(possibly a 1 element list) and process that way
			
			for parent_id in parent_id_base.split(","):
				
				if parent_id in trans_to_parent:
					parent_gene = trans_to_parent[parent_id]
				elif parent_id in genes:
					parent_gene = parent_id
					parent_id = "DUMMY_" + parent_id
					genes[parent_gene].add_trans(parent_id, Transcript(parent_id, parent_gene, "DUMMY"))
				else:
					e.write(' * PARENT_ERROR ' + line + "\n")
			
				#now parent_gene has the parent gene id, and parent_id is either a transcript or a gene
			
				if gff_fields[2] == "CDS":
					#add protein id for CDS features
					genes[parent_gene].update_trans(parent_id, "protid", elem_id)
					#also update biotype of parent -- CDS == protein_coding by definition
					genes[parent_gene].biotype = "protein_coding"
					
				genes[parent_gene].add_length(parent_id, elem_len, gff_fields[2])
			
#DONE PROCESSING			

if args.verbose:
	sys.stderr.write ("Done reading GFF file.\n")

##WRITE PICKLE FILE##

if pickle_file != "/dev/null":
	if args.verbose:
		sys.stderr.write ("Dumping gene dictionary to pickle file.\n")
	pickle.dump(genes, open(pickle_file, "wb"))

##WRITE TABLE OUTPUT###

if table_file != "/dev/null":
	if args.verbose:
		sys.stderr.write ("Writing table output.\n")
	
	#want to write a tab-deliminted file with the following columns:
	#gene id, gene acc, gene name, gene biotype, and the for each transcript
	#transcript id, transcript acc, transcript length, protein length, protein id, is longest protein?
	
	with open (table_file, 'w') as tab:
		tab.write("#" + "\t".join(["species", "gene_id", "gene_acc", "gene_name", "biotype", "trans_id", "trans_acc", "trans_len", "cds_len", "prot_id", "is_longest_prot"]) + "\n")
		for geneid in genes:
			curgene=genes[geneid]
			for key in curgene.transcripts:
				curtrans=curgene.get_trans(key)
				tab.write(species + "\t" + "\t".join([curgene.id, xstr(curgene.acc), xstr(curgene.name), xstr(curgene.biotype)]) + "\t") #gene info
				tab.write("\t".join([curtrans.id, xstr(curtrans.acc), xstr(curtrans.len), xstr(curtrans.cdslen), xstr(curtrans.protid)]) + "\t") #transcript info
				tab.write(curgene.is_longest_prot(curtrans.id) + "\n")
 
			
	
