# -*- coding: utf-8 -*-
"""
Created on Thu Oct 15 11:40:12 2015

@author: Phil
note that first strand will always be chicken.
"""

import csv
import re

#file_out = '/Users/Phil/Desktop/Chicken_messingaround.txt'
#fout = open(file_out, 'w')
with open('/Users/Phil/Desktop/fulGla.psl', 'rU') as handle: #opens bed in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader:
        name = strLine[0]
        strand = strLine[9]
        qName = strLine[10]
        qSize = int(strLine[11])
        qStart = int(strLine[12])
        qEnd = int(strLine[13])
        tName = strLine[14]
        tSize = int(strLine[15])
        tStart = int(strLine[16])
        tEnd = int(strLine[17])
        blockCount = int(strLine[18])
        blockSizes = strLine[19]
        qStarts = strLine[20]
        tStarts = strLine[21]
        qExonSize=qEnd-qStart
        tExonSize=tEnd-tStart        
        #print name+"_"+str(qExonSize)+"_"+str(tExonSize)
        if qExonSize != tExonSize:
            difExon = qExonSize-tExonSize            
            print name+"_"+str(difExon)
        #print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name,strand,qName,qSize,qStart,qEnd,tName,tSize,tStart,tEnd,blockCount,blockSizes,qStarts,tStarts)
        
        
        
        
        """Start = int(strLine[3])-1 #BED is 0 based and GFF is 1 based
        Stop = strLine[4]
        MetaD = strLine[8]
        Zero = 0
        Strand = strLine[6]
        IDsearch = re.search('Dbxref\=([a-zA-Z0-9\.\,\_\:\/\-]*)',MetaD)
        MetaData= IDsearch.group(1)+",Start:"+str(Start)+",Stop:"+str(Stop)+",Strand:"+Strand
        outStr = "%s\t%s\t%s\t%s\t%s\t%s" % (Scaffold,Start,Stop,MetaData,Zero,Strand)
        fout.write(outStr + '\n') #and write that string out
fout.close()"""