import csv
import re
import Bio
import os
import sys
from Bio import SeqIO

matches = set() #here we use a set because it only records unique items
#fout=open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/Incomplete_Files/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+os.getenv('FILE'), 'w')
with open('/n/home12/pgrayson/regal/PseudoSearch/Final/output/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+os.getenv('FILE'), 'rU') as handle: #opens maf in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        if len(strLine)>0: #if this line is greater than 0
            if strLine[0]=='a' and len(matches)!=0: #if it has an "a" at position 0
                lsttxt=[]
                for thing in matches: #for every item in matches
                    scaf_search = re.search('([a-zA-Z]*)\.([a-zA-Z0-9_.]*)',thing) #split that species.scaffold on the "."
                    genome = scaf_search.group(1) #call the species side genome #call the scaffold side scaff
                    lsttxt.append(genome)
                if len(set(lsttxt))!=len(matches):
                    fout=open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/Incomplete_Files/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+os.getenv('FILE'), 'w')
                    fout.write('duplicate')
                matches = set()
		break
            if strLine[0]=='s':
                matches.add(strLine[1])
fout.close()
