# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: Phil
"""

from Bio import SeqIO
from Bio.Seq import Seq

infile = '/Users/Phil/Desktop/fulGla_concat.txt'
for record in SeqIO.parse(infile,"fasta"):
    trec = record.seq.translate()
    strec = trec.split('*')
print trec
print strec
for rec in strec:
    print rec

                    
infile = '/Users/Phil/Desktop/test.txt'
for record in SeqIO.parse(infile,"fasta"):
    for strand, nuc in [(+1, record.seq)]:
        for pro in nuc.translate().split("*"):
            if len(pro) >= 1 and pro != len(pro) * "X":                
                print("%s - length %i, %s" % (pro[:30], len(pro), record.description))





