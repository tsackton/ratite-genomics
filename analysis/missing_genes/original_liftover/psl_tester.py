import csv
import os
import re
total = 0 #a bunch of counts that start at 0
perfect = 0
imperfect = 0
falseStart = 0
falseEnd = 0
falseBoth = 0

fout = open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/counts/'+os.getenv('FILE')+'.txt','w') #out files will be ${FILE}.txt
with open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/'+os.getenv('FILE'), 'rU') as handle:
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in and give names to bits that will be used
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
        qnamestart = re.search('Start\:([0-9]*)\,', name) #pull the chicken start position from name
        qnamestop = re.search('Stop\:([0-9]*)\,', name) #pull the chicken start position from name
        total += 1 #as you read a line, add 1 to total
        if int(qnamestart.group(1)) == qStart and int(qnamestop.group(1)) == qEnd: #if the chicken positions in name match those used in qStart and qEnd
            qExonSize=qEnd-qStart #check qExon size
            tExonSize=tEnd-tStart #check tExon size
            if qExonSize==tExonSize: #if both these exons are the same size
                perfect +=1 #this is a perfect sequence
            else:
                imperfect +=1 #otherwise, it is imperfect (i.e., it has indels)
        else: #if the chicken start and end don't jive, there are 3 more categories
            if int(qnamestart.group(1)) != qStart and int(qnamestop.group(1)) == qEnd: #missing start, but end is correct
                falseStart +=1
            elif int(qnamestop.group(1)) != qEnd and int(qnamestart.group(1)) == qStart: #missing end, but start is correct
                falseEnd +=1
            elif int(qnamestop.group(1)) != qEnd and int(qnamestart.group(1)) != qStart: #missing start and end
                falseBoth +=1
components = (int(perfect+imperfect+falseStart+falseEnd+falseBoth)) #add all the components up
fout.write("total = "+str(total)) #write the total
fout.write("\n"+"perfect = "+str(perfect)) #and each individual component
fout.write("\n"+"imperfect = "+str(imperfect))
fout.write("\n"+"falseStart = "+str(falseStart))
fout.write("\n"+"falseEnd = "+str(falseEnd))
fout.write("\n"+"falseBoth = "+str(falseBoth))
fout.write("\n"+"total is "+str(total)+" and components add to "+str(components)+"\n")
if int(total) != int(components): #if the total and the sum of components are not equal
    fout.write("\n"+"We've got a problem here"+"\n") #write this message
fout.close()
