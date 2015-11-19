# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:30:19 2015
This file reads through an NCBI gff and creates a new gff with only the lines that pertain to the CDS of protein coding genes.
@author: phil-grayson
"""

import csv
import sys
import os
import argparse

#better option processing
parser = argparse.ArgumentParser(description="Read a GFF file and export lines with a specified type.")
parser.add_argument('gff', help='GFF file to process', action='store')
parser.add_argument('-o', '--out', help='File to write to', action='store')
parser.add_argument('-t', '--type', help='GFF line type to get (default=CDS); for multiple fields use a comma-delimited list', default='CDS', action='store')
args=parser.parse_args()

matchTypes = args.type.split(",") #make a list of match types

#outfile
try:
	outputFile = args.out
except:
	outputFile = args.gff +"/out/" + args.gff.split('.')[0]+ ".gff"
	
fout = open(outputFile, 'w')
with open(args.gff, 'r') as handle: #opens gff in read mode mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        if len(strLine)==9: #only the lines in proper gff format (length of 9) should be examined
            if strLine[2] in matchTypes: #if a given line in file matches types in 3 ...
                fout.write('\t'.join(strLine) + '\n') #and write that string out
fout.close()
