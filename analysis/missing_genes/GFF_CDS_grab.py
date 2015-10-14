# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:30:19 2015
This file reads through an NCBI gff and creates a new gff with only the lines that pertain to the CDS of protein coding genes.
@author: phil-grayson
"""

import csv

file_out = '/Users/Phil/Desktop/PseudoSearch/Chicken_CDS.gff'
fout = open(file_out, 'w')
with open('/Users/Phil/Desktop/PseudoSearch/GCF_000002315.3_Gallus_gallus-4.0_genomic.gff', 'rU') as handle: #opens gff in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        if len(strLine)==9: #only the lines in proper gff format (length of 9) should be examined
            if strLine[2]=="CDS": #if a given line in file has CDS in column 3 ...
                strLine2 = '\t'.join(strLine) #turn that line (which is in list format) into a tab delimited string
                fout.write(strLine2 + '\n') #and write that string out
fout.close()