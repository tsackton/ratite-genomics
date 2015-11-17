# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 14:35:27 2015

NCBI formatted gff to bed converter.  Takes a file with each exon on separate line and produces bed with each exon on separate line.

@author: phil-grayson
"""

import csv
import re
import os
import sys

phile = str(sys.argv[1])

fout = open(os.getcwd()+"/out/"+phile.split('.')[0]+".bed", 'w')
with open(os.getcwd()+"/"+phile, 'rU') as handle: #opens bed in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        Scaffold = strLine[0]
        Start = int(strLine[3])-1 #BED is 0 based and GFF is 1 based
        Stop = strLine[4]
        MetaD = strLine[8]
        Zero = 0
        Strand = strLine[6]
        IDsearch = re.search('Dbxref\=([a-zA-Z0-9\.\,\_\:\/\-]*)',MetaD)
        try:
            MetaData = IDsearch.group(1)+",Start:"+str(Start)+",Stop:"+str(Stop)+",Strand:"+Strand
        except:
            MetaData = str(MetaD)+"Start:"+str(Start)+",Stop:"+str(Stop)+",Strand:"+Strand
        outStr = "%s\t%s\t%s\t%s\t%s\t%s" % (Scaffold,Start,Stop,MetaData,Zero,Strand)
        fout.write(outStr + '\n') #and write that string out
fout.close()