# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: Phil
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

infile = '/Users/Phil/Desktop/test.txt'
for record in SeqIO.parse(infile,"fasta"):
    ts_ed = [] #trans_split (empty dropped) 
    trans_split = record.seq.translate().split('*')
    print trans_split
    print len(trans_split)    
    if len(trans_split) == 2:
        trans_split.pop(1)
        print trans_split
    else:
        end = len(trans_split)
        if len(trans_split[end-1]) == 0:
            print trans_split

    for orf in trans_split:      
        norf = orf+star
        ts_ed.append(norf)
    print ts_ed 
    print trans_split       
        
        if len(orf) > 0:
            ts_ed.append(orf)
    end = len(ts_ed)
    star = "*"
    star = (Seq(star))
    nts_ed = []    
    for norf in ts_ed[0:end-1]:
        norf = norf+star
        nts_ed.append(norf)
    print nts_ed
    

infile = '/Users/Phil/Desktop/test.txt'
for record in SeqIO.parse(infile,"fasta"):
    ts_ed = [] #trans_split (empty dropped) 
    trans_split = record.seq.translate().split('*')
    for orf in trans_split:      
        if len(orf) > 0:
            ts_ed.append(orf)
    end = len(ts_ed)
    star = "*"
    star = SeqRecord(Seq(star))
    nts_ed = []    
    for norf in ts_ed[0:end-1]:
        norf = norf+star
        nts_ed.append(norf)
    print nts_ed






    
    for orf in ts_ed
    print max(mylist, key=len)

print trec
print strec
for rec in strec:
    print rec

