# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17, 2015

Convert the chicken liftover psl files into bed files.  

@author: phil-grayson
"""

import csv
import os
import sys

phile = str(sys.argv[1])

fout = open(os.getcwd()+"/bed/"+phile.split('.')[0]+"psl.bed", 'w')
with open(os.getcwd()+"/"+phile, 'rU') as handle: #opens bed in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        Scaffold = strLine[14]
        Start = int(strLine[16])
        Stop = strLine[17]
        MetaD = strLine[0]
        Zero = 0
        Strand = strLine[9][1:]
        outStr = "%s\t%s\t%s\t%s\t%s\t%s" % (Scaffold,Start,Stop,MetaD,Zero,Strand)
        fout.write(outStr + '\n') #and write that string out
fout.close()
