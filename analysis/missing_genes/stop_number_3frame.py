# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 17:10:11 2015

@author: phil-grayson
"""
from Bio import SeqIO
import re
import sys

phile = str(sys.argv[1])
nameo = phile.split("_")[1].split(".")[0]
fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/orfs/'+nameo+'_orfs.txt','w')
fout.write("ChickenSeqID"+"\t"+nameo+"\n")
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
    recordgrab = re.search('ChickenGeneID\:([0-9]*.*Genbank\:[A-Za-z0-9/./_]*)',record.description)
    fout.write(str(recordgrab.group(1))+"\t"+str(orfs)+"\n")    
fout.close()
