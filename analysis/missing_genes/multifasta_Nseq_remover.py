# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 14:03:25 2015

@author: phil-grayson
"""

import sys
import re
import csv
import os
from Bio import SeqIO

fout = open('/Users/Phil/Desktop/???','w')

infile = '/Users/Phil/Desktop/concat_???.fa'
for record in SeqIO.parse(infile,"fasta"):
    if str(record.seq) != len(record.seq) * 'N':
        fout.write(">"+str(record.description)+"\n"+str(record.seq)+"\n")
fout.close()