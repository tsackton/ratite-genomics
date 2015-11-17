# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:30:19 2015
This file reads through an NCBI gff and creates a new gff with only the lines that pertain to the CDS of protein coding genes.
@author: phil-grayson
"""

import csv
import sys
import os

phile = str(sys.argv[1])
fout = open(os.getcwd()+"/out/"+phile.split('.')[0]+".gff", 'w')

with open(os.getcwd()+"/"+phile, 'rU') as handle: #opens gff in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        if len(strLine)==9: #only the lines in proper gff format (length of 9) should be examined
            if strLine[2]=="CDS": #if a given line in file has CDS in column 3 ...
                fout.write('\t'.join(strLine) + '\n') #and write that string out
fout.close()
