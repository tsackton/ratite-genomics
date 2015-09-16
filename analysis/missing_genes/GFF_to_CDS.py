# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 13:30:19 2015
This file reads through an NCBI gff and creates a new gff with only the lines that pertain to the CDS of protein coding genes.
@author: phil-grayson
"""

import csv
import re

file_out = '/Users/Phil/Desktop/PseudoSearch/Chicken_CDS_Only.gff'
fout = open(file_out, 'w')
with open('/Users/Phil/Desktop/PseudoSearch/GCF_000002315.3_Gallus_gallus-4.0_genomic.gff', 'rU') as handle: #opens gff in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for lll in reader:
        if len(lll)==9: #only the lines in proper gff format (length of 9) should be examined
            pcsearch = re.search(r'GeneID\:(\d*)\;.*gene_biotype=protein_coding', lll[8]) #find lines with biotype=protein_coding and records as pcsearch
            if bool(pcsearch): #if protein coding biotype in line
                pcmatch = pcsearch.group(1) #call pcmatch the gene id of that gene
            if lll[2]=="CDS" and str(pcmatch) in lll[8]: #if a given line in file has CDS in column 3 and gene id matches in column 9...
                llll = '\t'.join(lll) #turn that line (which is in list format) into a tab delimited string
                fout.write(llll + '\n') #and write that string out
fout.close()

