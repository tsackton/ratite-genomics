# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: phil-grayson
"""
from Bio import SeqIO
from Bio.Seq import Seq

star = "*"
star = (Seq(star))
ts_ed = []
fout = open('/Users/Phil/Desktop/testoutgG_orfs.txt','w')
fout2 = open('/Users/Phil/Desktop/testoutgG_longest.txt','w')
fout3 = open('/Users/Phil/Desktop/testoutgG_trans.txt','w')
infile = '/Users/Phil/Desktop/concat_galGal.fa'
for record in SeqIO.parse(infile,"fasta"):
    trans = record.seq.translate()    
    trans_split = trans.split('*') #translate the sequence and split on stops  
    orfs = len(trans_split)
    if len(trans_split[-1]) != 0 and len(trans_split[0]) != 0: #if there was not a stop at the beginning or end of the translated protein
        rf = trans_split[-1] #call the last reading frame rf
    if len(trans_split[0]) == 0: #if there are any stop codons at the beginning of the sequences
        try: #necessary because an entire sequence of stop codons would throw an index error
            while len(trans_split[0]) == 0: #while there are stop codons at the beginning of the sequence
                trans_split.pop(0) #pop them off
            rf = trans_split[-1] #once the stops are gone, call the last reading frame as rf
        except:
            rf = "*" #if the entire sequence was stop codons, call the last reading frame "*"
    try: #as above, necessary in case the whole sequence was stop codons
        if len(trans_split[-1]) == 0: #if you have stop codons at the end of the sequence
            while len(trans_split[-1]) == 0: #while there are stop codons at the end of the sequence
                trans_split.pop(-1) #pop them off
            rf = trans_split[-1]+star #once the stops are gone, add one back to the final piece of sequence
    except:
        None #do nothing since rf is already "*" from above
    for orf in trans_split[0:-1]: #for every reading frame for a given transcript (except the final one)
        orf = orf+star #add "*" back, since we've split it off
        ts_ed.append(orf) #to a new list, add these pieces in order
    ts_ed.append(rf) #and once that is complete, add on the final reading frame (either "*", seq with "*", or seq without "*" based on biology)
    longest = max(ts_ed, key=len) #call the longest reading frame for a given transcript longest
    ts_ed=[]
    fout.write(str(record.description)+"\t"+str(orfs)+"\n")    
    fout2.write(">"+str(record.description)+"\n"+str(longest)+"\n") #and write this out with appropriate records
    fout3.write(">"+str(record.description)+"\n"+str(trans)+"\n")
fout.close()
fout2.close()
fout3.close()