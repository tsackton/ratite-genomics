# -*- coding: utf-8 -*-
"""
Created on Mon Feb  8 17:31:41 2016

@author: Phil
"""

from Bio import SeqIO
import csv

fout = open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/emu_transcripts_test.txt', 'w') #outfile
with open('/Users/Phil/Dropbox/Doctorate/Tabin_Limbs/emu_primer_work/emu_first_run_with_gene_names.txt', 'rU') as handle: #infile
    reader=csv.reader(handle,delimiter='\t')
    for strLine in reader:
        for record in SeqIO.parse('/Users/Phil/Desktop/droNov.genome.all.maker.transcripts.fasta','fasta'):
            if strLine[0] in record.description:               
                fout.write(">"+str(record.description).split("-")[0]+'_'+strLine[1]+'\n'+str(record.seq)+"\n")
fout.close()