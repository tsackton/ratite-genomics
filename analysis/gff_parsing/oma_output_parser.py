# takes input like: python hogget_working1.py [directory of HOGs] [protein list]

import sys, re, os, string, exceptions, itertools, pprint, os.path
from collections import OrderedDict
from itertools import groupby
from operator import itemgetter

# species name
# species6 = sys.argv[1]
# input directory
input = sys.argv[1]
# input source list of proteins
proteinlist = sys.argv[2]

species6 = proteinlist[:6]

HOGassigned = species6 + '.assigned'
noHOG = species6 + '.unassigned'

hogs = os.listdir(input)
hogassigned = []

class Prot(object):
	def __init__(self, geneID, proteinID, star):
		self.geneID = geneID
		self.proteinID = proteinID
		self.star = star

	def __hash__(self):
		return hash((self.geneID, self.proteinID, self.star))

	def __eq__(self, other):
		return (self.geneID, self.proteinID, self.star) == (other.geneID, other.proteinID, other.star)

for hog in hogs:
	path = input + '/' + hog
	with open (path) as f:
		lines = f.read().splitlines()
		for line in lines:
			if species6 in line:
				if not re.findall('([a-zA-Z]+P_\d+).', line):
					print 'PROT ' + line
					continue
				else:
					protein = re.findall('([a-zA-Z]+P_\d+).', line)[0]
                
				if not re.findall('\|(\d+)\|', line):
					print 'GENE ' + line
					geneID = ''
				else:
					geneID = re.findall('\|(\d+)\|', line)[0]
                   
				hogassigned.append(Prot(geneID, [protein, hog], ''))

filetowrite = []
proteinsWithHOG = []
listofproteins = []
listofgenes = []

with open (proteinlist) as f:
	for line in f:
		if not line.startswith(species6):
			continue
		if not re.findall('protein_id=(\S+_\d+).', line):
			print 'PROT 1 ' + line
			continue
		else:
			proteinToMatch = re.findall('protein_id=(\S+_\d+).', line)[0]
			geneToMatch = re.findall('GeneID:(\d+)', line)[0]
			if not '*' in line:
				thingToMatch = Prot(geneToMatch, proteinToMatch, '')
			else:
				thingToMatch = Prot(geneToMatch, proteinToMatch, '*')
		if thingToMatch.geneID in listofgenes:
#			print thingToMatch.geneID
			continue
		else:
			listofgenes.append(thingToMatch.geneID)
			listofproteins.append(thingToMatch)
		for hog in hogassigned:
			if hog.proteinID[0] == thingToMatch.proteinID:
				proteinsWithHOG.append(thingToMatch)
				line = line.replace("\n","")
				line += '\t' + hog.proteinID[1]
				filetowrite.append(line)
      #          print hog.geneID
       #         print thingToMatch.geneID
			elif hog.geneID == thingToMatch.geneID:
				proteinsWithHOG.append(thingToMatch)
				line = line.replace("\n","")
				line += '\t' + hog.proteinID[1] + ' *'
				filetowrite.append(line)

unassigned = list(set(listofproteins)-set(proteinsWithHOG))

with open(noHOG, 'w') as f:
	for protein in unassigned:
		f.write(protein.geneID + '\t' + protein.proteinID + '\t' + protein.star)
		f.write('\n')

with open(HOGassigned, 'w') as f:
	for line in filetowrite:
		f.write(line)
		f.write('\n')
