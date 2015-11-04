# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: phil-grayson
"""
from Bio import SeqIO
import sys

phile = str(sys.argv[1])
nameo = phile.split("_")[1].split(".")[0]
fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/stopnums/step1/'+nameo+'_into.txt','w')
for record in SeqIO.parse(phile,"fasta"):
    trans = record.seq.translate() #frame 1
    trans2 = record.seq[1:].translate() #frame 2
    trans3 = record.seq[2:].translate() #frame 3 
    trans_split = trans.split('*') #translate the sequence and split on stops 
    trans2_split = trans2.split('*')
    trans3_split = trans3.split('*')
    orfs = len(trans_split)-1
    orfs2 = len(trans2_split)-1
    orfs3 = len(trans3_split)-1
    frames = [orfs,orfs2,orfs3]
    orfs4 = min(frames) #best of 3 frames
    fout.write(record.description+"\t"+str(orfs)+"\t"+str(orfs2)+"\t"+str(orfs3)+"\t"+str(orfs4)+"\n")    
fout.close()
