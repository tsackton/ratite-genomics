# -*- coding: utf-8 -*-
"""
Created on Fri Oct 30 12:17:41 2015

@author: phil-grayson
"""
import sys
import re
import csv

phile = str(sys.argv[1])
nameo = phile.split("_")[0].split(".")[0]

total = 0
overall_2 = 0
overall_2plus = 0
minusminus_2 = 0
minusminus_2plus = 0
minusplus_2 = 0
minusplus_2plus = 0
plusminus_2 = 0
plusminus_2plus = 0
plusplus_2 = 0
plusplus_2plus = 0
none_2 = 0
none_2plus = 0
fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/stopnums/'+nameo+'best.txt','w')

with open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/stopnums/step1/'+phile, 'rU') as handle:
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #line by line, calls the first column name
        total += 1
        des = strLine[0]
        orfs = strLine[4] #1 gives frame 1, 4 gives best of 3
        strand = re.search('pslStrand\:([\-\+None]*)\,',des)
        if int(orfs) == 1 or int(orfs) == 0:
            overall_2 += 1
            if str(strand.group(1)) == '--':
                minusminus_2 += 1
            elif str(strand.group(1)) == '-+':
                minusplus_2 += 1
            elif str(strand.group(1)) == '+-':
                plusminus_2 += 1
            elif str(strand.group(1)) == '++':
                plusplus_2 += 1
            elif str(strand.group(1)) == 'None':
                none_2 += 1
            else:
                print "we have a problem with 2's"+strLine
        else:
            overall_2plus += 1
            if str(strand.group(1)) == '--':
                minusminus_2plus += 1
            elif str(strand.group(1)) == '-+':
                minusplus_2plus += 1
            elif str(strand.group(1)) == '+-':
                plusminus_2plus += 1
            elif str(strand.group(1)) == '++':
                plusplus_2plus += 1
            elif str(strand.group(1)) == 'None':
                none_2plus += 1
            else:
                print "we have a problem with 2plus's"+strLine

components = (int(minusminus_2+minusminus_2plus+minusplus_2+minusplus_2plus+plusminus_2+plusminus_2plus+plusplus_2+plusplus_2plus+none_2+none_2plus))
fout.write("total = "+str(total)) #write the total
fout.write("\n"+"overall_2 = "+str(overall_2)) #and each individual component
fout.write("\n"+"overall_2plus = "+str(overall_2plus))
fout.write("\n"+"minusminus_2 = "+str(minusminus_2))
fout.write("\n"+"minusminus_2plus = "+str(minusminus_2plus))
fout.write("\n"+"minusplus_2 = "+str(minusplus_2))
fout.write("\n"+"minusplus_2plus = "+str(minusplus_2plus))
fout.write("\n"+"plusminus_2 = "+str(plusminus_2))
fout.write("\n"+"plusminus_2plus = "+str(plusminus_2plus))
fout.write("\n"+"plusplus_2 = "+str(plusplus_2))
fout.write("\n"+"plusplus_2plus = "+str(plusplus_2plus))
fout.write("\n"+"none_2 = "+str(none_2))
fout.write("\n"+"none_2plus = "+str(none_2plus))
fout.write("\n"+"total is "+str(total)+" and components add to "+str(components)+"\n")
fout.write("\n"+"overall, "+str("{0:.2f}".format((float(overall_2plus)/total)*100))+" percent of the transcripts have more than 2 orfs"+"\n")
fout.close()
