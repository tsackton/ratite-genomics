import csv
import os
import re
import sys
IDic = {}
single = []
duplicate = []

fout_singles = open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/singles/singles_'+os.getenv('FILE'),'w') #out files will be ${FILE}.txt
fout_duplicates = open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/duplicates/duplicates_'+os.getenv('FILE'),'w') 
with open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/'+os.getenv('FILE'), 'rU') as handle: #opens psl in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in and give names to bits that will be used
        name = strLine[0]
        namegrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,(.*\:[\+\-])',name)
        IDic.setdefault(str(namegrab.group(1)),[]).append(str(namegrab.group(2)))
    for key in IDic:
        if len(set(IDic[key]))==len(IDic[key]):
            single.append(key)
        else:
            duplicate.append(key)
with open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/'+os.getenv('FILE'), 'rU') as handle:
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in and give names to bits that will be used
        name = strLine[0]        
        namegrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,(.*\:[\+\-])',name)       
        if namegrab.group(1) in duplicate:
            strLine2 = '\t'.join(strLine) #turn that line (which is in list format) into a tab delimited string
            fout_duplicates.write(strLine2 + '\n')             
        else:
            strLine2 = '\t'.join(strLine) #turn that line (which is in list format) into a tab delimited string
            fout_singles.write(strLine2 + '\n')                        
    fout_duplicates.close()
    fout_singles.close()
