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
fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/orfs/'+nameo+'_stops.txt','w')
fout.write("ChickenSeqID"+"\t"+nameo+"_frame1"+"\t"+nameo+"_frame2"+"\t"+nameo+"_frame3"+"\t"+nameo+"_best"+"\n")
fout2 = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/translation/'+nameo+'_trans3frames.txt','w')
fout3 = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/translation/'+nameo+'_transbestframe.txt','w')
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
    frames = [orfs,orfs2,orfs3]
    orfs4 = min(frames)
    fout2.write(">frame1 "+str(record.description)+"\n"+str(trans)+"\n"+">frame2 "+str(record.description)+"\n"+str(trans2)+"\n"+">frame3 "+str(record.description)+"\n"+str(trans3)+"\n")
    if orfs4 == orfs:
	fout3.write(">bframe1_"+nameo+"_"+str(record.description)+"\n"+str(trans)+"\n")
    elif orfs4 == orfs2 == orfs3:
        fout3.write(">bframe2_"+nameo+"_"+str(record.description)+"\n"+str(trans2)+"\n"+">bframe3_"+nameo+"_"+str(record.description)+"\n"+str(trans3)+"\n")
    elif orfs4 == orfs2:
        fout3.write(">bframe2_"+nameo+"_"+str(record.description)+"\n"+str(trans2)+"\n")
    elif orfs4 == orfs3:
        fout3.write(">bframe3_"+nameo+"_"+str(record.description)+"\n"+str(trans3)+"\n")
    if str(record.seq) == len(record.seq) * 'N':
        orfs4 = -1
    fout.write(str(recordgrab.group(1))+"\t"+str(orfs)+"\t"+str(orfs2)+"\t"+str(orfs3)+"\t"+str(orfs4)+"\n")
fout.close()
fout2.close()
fout3.close()
