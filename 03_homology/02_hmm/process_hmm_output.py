#!/usr/bin/env python

import os, sys, argparse, string, re
import networkx as nx

SINGLE_INDEX=1
PGRAPH=nx.DiGraph()

def process_hog(inputdir, fasta_ext):
	protdict={}
	for hog in os.listdir(inputdir):
		if hog.endswith(fasta_ext):
			hog_to_open = inputdir + '/' + hog
		
			#hog name trims off the extension
			hog_name=hog.strip(("." + fasta_ext))

			#all we care about are parsing deflines, so don't bother with the overhead of biopython

			with open (hog_to_open, 'r') as f:
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
					
						protdict[cur_prot]=hog_name
	return(protdict)

def process_info(inputfile):
	protinfo={}
	with open (inputfile, "r") as ps:
		for line in ps:
			line.rstrip("\n")
			if line.startswith("#"):
				continue
		
			fields=line.split()
			try:
				fields[9] = re.sub('\.\d+$', '', fields[9])
			except:
				continue # apply your error handling
		
			found=fields[9]		
			protinfo[found]={}
			protinfo[found]['sp']=fields[0]
			protinfo[found]['gene']=fields[2]
	return protinfo

def process_hmmsearch(hmminput):
	hmmsearch={}
	with open (hmminput, 'r') as hmm:
		for line in hmm:
			hmmFields=line.split()
			#fields[0] = protein, fields[2] = hog, fields[4]  = eval, fields[5] = bitscore
		
			if hmmFields[0] in hmmsearch:
				#already seen this protein, so now need to check if the evalue is at minimum
				cur_max=float(hmmsearch[hmmFields[0]]['score'])
				if float(hmmFields[5]) > cur_max:
					hmmsearch[hmmFields[0]]['hogs']=[hmmFields[2]]
					hmmsearch[hmmFields[0]]['score']=hmmFields[5]
				elif hmmFields[5] == cur_max:
					hmmsearch[hmmFields[0]]['hogs'].append(hmmFields[2])
				else:
					continue
			else:
				#have not seen this protein, need to create dict entry
				hmmsearch[hmmFields[0]]={}
				hmmsearch[hmmFields[0]]['score']=hmmFields[5]
				hmmsearch[hmmFields[0]]['hogs']=[]
				hmmsearch[hmmFields[0]]['hogs'].append(hmmFields[2])
	return hmmsearch
	
def print_hog_matrix(inmat, sp, outfile):
	allsp=sp
	hogmat=inmat
	
	with open (outfile, 'w') as out:
		header = "hog\t" + "\t".join(sorted(allsp))
		print(header, file=out)

		for hogid in hogmat:
			print(hogid, end="", file=out)
			for sp in sorted(allsp):
				try:
					spct = hogmat[hogid][sp]
				except KeyError:
					spct = 0
				
				print("\t" + str(spct), end="", file=out)
			print("\n", end="", file=out)

## PROCESS OPTIONS ##
parser = argparse.ArgumentParser(description="Read a directory of OMA output and hmm search output and new HOG matrix.")
parser.add_argument('hmm', help='HMM search output', action='store')
parser.add_argument('hog', help='HOG Directory', action='store')
parser.add_argument('info', help="Info file", action='store')
parser.add_argument('-f', '--fasta_ext', help="Extension of the fasta files in the HOG dir, defaults to fa", action="store", default="fa")
parser.add_argument('-p', '--protein', help='List of all proteins', action='store', default='hmm_prot_list.txt')
parser.add_argument('-s', '--species', help='List of all species', action='store', default='sp_list.txt')
parser.add_argument('-g', '--graph', help="Pickled graph file to start from", action='store')

args=parser.parse_args()

if not args.graph:
	#get dictionary mapping proteins to HOGs, singletons are not included
	prot_to_hog=process_hog(args.hog, args.fasta_ext)

	#make hog to count data
	hog_ct={}
	for prot in prot_to_hog:
		hog=prot_to_hog[prot]
		try:
			hog_ct[hog] += 1
		except KeyError:
			hog_ct[hog] = 1

	#process info
	protinfo=process_info(args.info)

	#make species set
	allsp=set(open(args.species).read().split())

	#proces hmm search file
	hmmsearch=process_hmmsearch(args.hmm)

	## GET PROTEIN LIST ##
	with open (args.protein, 'r') as prot:
		protlist=[line.rstrip('\n') for line in prot]

	## MAIN WORK: MAKE GRAPH

	for protein in protlist:
		if not protein in prot_to_hog:
			new_hog_id = "HOGSingle" + str(SINGLE_INDEX)
			prot_to_hog[protein] = new_hog_id
			hog_ct[new_hog_id] = 1
			SINGLE_INDEX += 1
		
		init_hog=prot_to_hog[protein]

		if not init_hog in PGRAPH:
			PGRAPH.add_node(init_hog)
		
		try:
			hmm_hogs=hmmsearch[protein]['hogs']
		except:
			hmm_hogs=[]
	
		try:
			PGRAPH.node[init_hog]['prot'].append(protein)
		except:
			PGRAPH.node[init_hog]['prot'] = [ protein ]
	
		#hmm_hogs is a (possibly empty) list with all the hmm hits for that protein
		#for each element in the hmm_hog, add an edge from init_hog to hmm_hog, updated edge count if necessary
		for new_hog in hmm_hogs:
			if new_hog == init_hog:
				continue
			
			try:
				PGRAPH[init_hog][new_hog]['count'] += 1
			except:
				PGRAPH.add_edge(init_hog, new_hog, count=1)
		
			PGRAPH[init_hog][new_hog]['prop_init'] = PGRAPH[init_hog][new_hog]['count'] / hog_ct[init_hog]
			PGRAPH[init_hog][new_hog]['prop_new'] = PGRAPH[init_hog][new_hog]['count'] / hog_ct[new_hog]
	
	#write graph in pickle format
	nx.write_gpickle(PGRAPH, "protein_graph.gpickle")

	#process graph
	hogmat={}
	
	for hog in PGRAPH.nodes():
		for prot in PGRAPH.node[hog]['prot']:
			try:
				this_spe = protinfo[prot]['sp']
				this_gene = protinfo[prot]['gene']
			except KeyError:
				if prot == "YP_007890942":
					this_spe = "pygAde"
				else:
					print(hog + " failed with prot " + prot, file=sys.stderr)
					continue
		
			try:
				hogmat[hog][this_spe] += 1
			except KeyError:
				if not hog in hogmat:
					hogmat[hog]={}
				hogmat[hog][this_spe]=1

	print_hog_matrix(hogmat, allsp, "original_hog_matrix.txt")
else:
	PGRAPH=nx.read_gpickle(args.graph)

#update
PNEW=nx.DiGraph(PGRAPH)
for edge in PNEW.edges(data="count"):
	size_test=min(len(PNEW.node[edge[0]]['prot']),len(PNEW.node[edge[1]]['prot']))
	if size_test > 20 and float(edge[2])/float(size_test) > 0.2:
		threshold=0
	elif size_test <= 20 and float(edge[2])/float(size_test) > 0.05:
		threshold=0
	else:
		threshold=1
	
	if threshold:
		PNEW.remove_edge(*edge[:2])

#merge weakly connected components
PSUB=sorted(nx.weakly_connected_component_subgraphs(PNEW), key=len, reverse=True)
PUPDATE=nx.Graph()
INDEX=1
for G in PSUB:
	new_name = "HOG2_" + str(INDEX)
	PUPDATE.add_node(new_name)
	INDEX += 1
	for node in G:
		try:
			PUPDATE.node[new_name]['prot'].append(G.node[node]['prot'])
		except:
			PUPDATE.node[new_name]['prot'] = [G.node[node]['prot']]
			
hogmatnew={}

with open ("new_hog_list.txt", 'w') as nhog:
	prot_by_node=nx.get_node_attributes(PUPDATE,'prot')
	for hog in prot_by_node:
		for prot_list in prot_by_node[hog]:
			for prot in prot_list:
				try:
					this_spe = protinfo[prot]['sp']
					this_gene = protinfo[prot]['gene']
				except KeyError:
					if prot == "YP_007890942":
						this_spe = "pygAde"
						this_gene = "NA"
					else:
						print(hog + " failed with prot " + prot, file=sys.stderr)
						continue
		
				try:
					hogmatnew[hog][this_spe] += 1
				except KeyError:
					if not hog in hogmatnew:
						hogmatnew[hog]={}
					hogmatnew[hog][this_spe]=1
				
				print(hog + "\t" + prot+ "\t" + this_gene + "\t" + this_spe, file=nhog)
	
print_hog_matrix(hogmatnew, allsp, "updated_hog_matrix.txt")

