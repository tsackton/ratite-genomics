import csv
import re
import Bio
import os
import sys
from Bio import SeqIO
matches = set() #here we use a set because it only records unique items
fout2 = open('/n/regal/edwards_lab/phil/PseudoSearch/genomes/raw/lists/'+'list'+os.getenv('SLURM_ARRAY_TASK_ID')+'.txt','w')

with open('/n/home12/pgrayson/regal/PseudoSearch/Final/output/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+os.getenv('FILE'), 'rU') as handle: #opens maf in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        if len(strLine)>0: #if this line is greater than 0        
            if strLine[0]=='s': #if it has an "s" at position 0
                matches.add(strLine[1]) #add the species.scaffold to the matches set
lsttxt=[]
for thing in matches: #for every item in matches
    scaf_search = re.search('([a-zA-Z]*)\.([a-zA-Z0-9_.]*)',thing) #split that species.scaffold on the "."
    genome = scaf_search.group(1) #call the species side genome
    scaff = scaf_search.group(2) #call the scaffold side scaff    
    lsttxt.append(genome+scaff+'.fa') #append genome+scaff with .fa to the list
    fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/temp/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+genome+scaff+'.fa', 'w') #this creates a unique file for each genome and scaffold
    file=genome+".fa" #this allows SeqIO to call the file from our current directory
    z = SeqIO.index_db(file+".idx", file, "fasta") #this calls the file and the already generated index file
    fout.write(z.get_raw(scaff)) #and writes the scaffold to fout
    fout.close()

for l in lsttxt:
    fout2.write('../temp/'+os.getenv('SLURM_ARRAY_TASK_ID')+'/'+l+',') #here we create a txt doc to feed to the m2fStitcher for the seq option with every genome and scaffold we've just loaded to the temp file
fout2.close()
