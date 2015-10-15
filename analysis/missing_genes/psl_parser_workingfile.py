# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:40:12 2015

@author: Phil
note that first strand will always be chicken.
"""

import csv
from Bio import SeqIO
import os

fout = open('/Users/Phil/Desktop/Chicken_messingaround.txt','w')
os.chdir('/Users/Phil/Desktop/')
with open('/Users/Phil/Desktop/fulGla.psl', 'rU') as handle: #opens bed in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        name = strLine[0]
        strand = strLine[9]
        qName = strLine[10]
        qSize = int(strLine[11])
        qStart = int(strLine[12])
        qEnd = int(strLine[13])
        tName = strLine[14]
        tSize = int(strLine[15])
        tStart = int(strLine[16])
        tEnd = int(strLine[17])
        blockCount = int(strLine[18])
        blockSizes = strLine[19]
        qStarts = strLine[20]
        tStarts = strLine[21]
        if strand == "+-":
            rc = True
        elif strand == "--":
            rc = True
        else:
            rc = False
        qExonSize=qEnd-qStart
        tExonSize=tEnd-tStart
        if qExonSize==tExonSize:
            #genome=os.getenv('1')#link to $1 in batch script
            #genome=str(file)+".fa"
            genome='fulGla.fa'
            z = SeqIO.index_db(genome+".idx", genome, "fasta")
            z = z.get(tName)[tStart:tEnd] #returns the target exon from the target scaffold as a SeqRecord 
            fasta = z.format("fasta") #returns the SeqRecord in fasta format
            if rc:
                zrc=(z.reverse_complement(id="rc_",description=z.description))
                rcfasta=zrc.format("fasta")
                fout.write(rcfasta)
            else:
                fout.write(fasta)
    fout.close()