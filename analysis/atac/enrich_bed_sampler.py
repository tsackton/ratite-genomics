# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 16:06:55 2017

argv annotation:

1 - number of replicates
2 - name of outfile
3 - name of peak file
4 - name of file to count lines from (accel file)
5 - name of file to sample from (all CNEE file or accel file)
6 - name of new outdir

@author: Phil
"""
total_overlap = [] #creates empty list that will be populated with lists of [total # bases sampled,total # bases overlapped in bed intersect]

import random
import subprocess
import sys

lines = sum(1 for line in open(sys.argv[4])) #sum number of lines in accelerated CNEE bed file

FileList = [] #creates empty list that will be populated with the lines read from the file we are sampling from (all CNEE or accel file)
with open(sys.argv[5]) as file:
    for line in file:
        line = line.strip() #strip newline
        FileList.append(line) 

subprocess.check_output("mkdir -p pyout_"+sys.argv[6],shell=True) #make unique directory for temp and outfiles

for i in range(int(sys.argv[1])): #loops for n times (e.g., 5000)
    boots = [] #empty list that will be populated with file lines sampled randomly from FileList
    for element in range(lines): #for each line in the accel file, sample a random line (with replacement) from FileList 
        boots.append(FileList[(random.randint(0,lines-1))])
    baseTot = 0 #start a count of total bases
    for boot in boots:
        #print baseTot
        start = int(boot.split("\t")[1]) #start coordinate of CNEE
        stop = int(boot.split("\t")[2]) #end coordinate of CNEE
        baseTot += (stop-start) #add the difference to baseTot as running tally
    
    
    prefx = str(sys.argv[3]).split("/")[-1] #used because argv[3] was called using ../gat/$FILE originally - just strips the $FILE off by itself
    fileout = ("pyout_"+sys.argv[6]+"/"+prefx+"_"+str(sys.argv[4])+"_"+str(sys.argv[5]+".txt")) #make a unique temp file name
    fout = open(fileout, 'w')
    fout.write('\n'.join(boots)) #write out each line from boots list
    fout.close()


    output = subprocess.check_output("bedtools intersect -a "+fileout+" -b "+sys.argv[3]+" -wo",shell=True) #call to bedtools intersect and return the output as output

    
    overlap = 0 #start an overlap count at 0
    outList = output.split("\n")[0:-1] #split output on new lines and remove the last one (empty) 
    for intersect in outList:
        overlap += int(intersect.split("\t")[-1]) #for every intersect in outList add the number of overlapping bases
    
    total_overlap.append([str(baseTot),str(overlap)]) #add the baseTot (total possible bases) and corresponding overlap to total_overlap list as new list of 2 elements

fileout = ("pyout_"+sys.argv[6]+'/'+sys.argv[2]) #new fileout (this is the final file)
fout = open(fileout, 'w')
for listel in total_overlap: #for every listel (2 element list) in total_overlap, write the elements separated by a tab and followed by a new line
    fout.write("\t".join(listel) + "\n")
        
fout.close()
