# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 14:23:31 2015

@author: Phil
"""

from Bio import SeqIO

fout = open('/Users/Phil/Desktop/3orfmesUni.txt','w')
infile = '/Users/Phil/Desktop/concat_mesUni.fa'
for record in SeqIO.parse(infile,"fasta"):
    trans = record.seq.translate()    
    trans_split = trans.split('*') #translate the sequence and split on stops  
    orfs = len(trans_split)
    if orfs == 3:  
        fout.write(">"+str(record.description)+"\n"+str(trans)+"\n")
fout.close()
