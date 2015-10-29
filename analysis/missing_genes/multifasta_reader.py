# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: phil-grayson
"""
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
star = "*"
star = (Seq(star))
ts_ed = []
nts_ed = []

infile = '/Users/Phil/Desktop/test.txt'
for record in SeqIO.parse(infile,"fasta"):
    ts_ed = [] #trans_split (empty dropped) 
    trans_split = record.seq.translate().split('*')      
    #print trans_split    
    if len(trans_split[-1]) == 0 and len(trans_split[0]) == 0:
        trans_split.pop(-1)
        trans_split.pop(0)
        rf = trans_split[-1]+star
    elif len(trans_split[-1]) == 0 and len(trans_split[0]) != 0:
        trans_split.pop(1)
        rf = trans_split[-1]+star
    elif len(trans_split[-1]) != 0 and len(trans_split[0]) == 0:
        trans_split.pop(0)           
        rf = trans_split[-1]
    else:
        rf = trans_split[-1]
    for orf in trans_split[0:-1]: 
        orf = orf+star
        nts_ed.append(orf)
    nts_ed.append(rf)
    print nts_ed
    nts_ed=[]
    



    for orf in ts_ed
    print max(mylist, key=len)


