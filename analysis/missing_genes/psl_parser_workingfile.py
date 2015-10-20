import csv
from Bio import SeqIO
import os
import re
total = 0 #a bunch of counts that start at 0
perfect = 0
imperfect = 0
falseStart = 0
falseEnd = 0
falseBoth = 0
x = [] #empty list for matches
y = [] #empty list to start collecting bed lines in  

#fout = open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/counts/'+os.getenv('FILE')+'.text','w') #out files will be ${FILE}.txt
#with open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/'+os.getenv('FILE'), 'rU') as handle:
fout = open('/Users/Phil/Desktop/messingaroundfulGla.txt','w')
os.chdir('/Users/Phil/Desktop/')
with open('/Users/Phil/Desktop/10000fulGla.psl', 'rU') as handle: #opens bed in universal mode
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
        qnamestop = re.search('Stop\:([0-9]*)\,', name) #pull the chicken end position from name
        if len(x)==0:  #if the program has just started, there is nothing in x      
            namegrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,(.*\:[\+\-])',name)                        
            genegrab = namegrab.group(1) #everything in re brackets becomes genegrab
            x = [genegrab] #and now x is genegrab too
            if x in name:
                y.append(strLine)
        if x[0] == str(genegrab):
            if x in name:
                y.append(strLine)
            if x not in name:
                #start doing something.
#                if str(pcmatch) not in strLine[3]: #if pcmatch is not in column 4 of a line (i.e., we have all exons to a transcript)
#                    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #write a file
#                    fout= open(fileout, 'w')
#                    while len(y)>0: #since y is a list of lists (one list for each exon in the transcript) we pop it until it is len==0
#                        y_item = (y.pop(0)) #pop the first list (bed line) off
#                        y_item_tab = '\t'.join(y_item) #tab seperate it
#                        fout.write(y_item_tab + '\n') #and write it to the file
#                    fout.close() #close that file -  we're done with it
#                    y = [] #start y over
#                    i += 1 #add 1 to the file counter
#                    y.append(strLine) #add the first line of the new transcript to y
#                    pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #same as above           
                    pcmatch = pcsearch.group(1)
                    x = [pcmatch]
            else: #same as above except
                if x in name:
                    y.append(strLine)
                if x not in name:
                #start doing something.
                
               
#                    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #write a file
#                    fout= open(fileout, 'w')
#                    while len(y)>0: #since y is a list of lists (one list for each exon in the transcript) we pop it until it is len==0
#                        y_item = (y.pop(0)) #pop the first list (bed line) off
#                        y_item_tab = '\t'.join(y_item) #tab seperate it
#                        fout.write(y_item_tab + '\n') #and write it to the file
#                    fout.close() #close that file.  we're done with it
#                    y = [] #start y over
#                    i += 1 #add 1 to the file counter
#                    y.append(strLine) #add the first line of the new transcript to y
#                    pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #same as above           
#                    pcmatch = pcsearch.group(1)
#                    x = [pcmatch]
    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #here because after every line has been read, we haven't popped y into a file
    fout= open(fileout, 'w')
    while len(y)>0: 
        y_item = (y.pop(0)) 
        y_item_tab = '\t'.join(y_item)
        fout.write(y_item_tab + '\n') 
    fout.close()       
        
        
#VERY ROUGH ABOVE.  This will need to be incorporated.       
        
        
            
                if str(pcmatch) == strLine[3]: #if the current line has pcmatch in column 4
                    y.append(strLine) #append that line to list y                 
        
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
            #else: #missing start and end
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

#            if strand == "+-":
#                rc = True
#            elif strand == "--":
#                rc = True
#            else:
#                rc = False
#            qExonSize=qEnd-qStart
#            tExonSize=tEnd-tStart
#            if qExonSize==tExonSize:
#                #genome=os.getenv('1')#link to $1 in batch script
#                #genome=str(file)+".fa"
#                genome='fulGla.fa'
#                z = SeqIO.index_db(genome+".idx", genome, "fasta")
#                z = z.get(tName)[tStart:tEnd] #returns the target exon from the target scaffold as a SeqRecord 
#                fasta = z.format("fasta") #returns the SeqRecord in fasta format
#                if rc:
#                    zrc=(z.reverse_complement(id="rc_",description=z.description))
#                    rcfasta=zrc.format("fasta")
#                    fout.write(rcfasta)
#                else:
#                    fout.write(fasta)
#    fout.close()
            
        
#        else:
#            #here is where we write the n's in.
#        if int(qnamestart.group(1)) != qStart:
#            print "Start is off for "+str(strLine)
#        if int(qnamestop.group(1)) != qEnd:
#            print "End is off for "+str(strLine)
#        
#        if strand == "+-":
#            rc = True
#        elif strand == "--":
#            rc = True
#        else:
#            rc = False
#        qExonSize=qEnd-qStart
#        tExonSize=tEnd-tStart
#        if qExonSize==tExonSize:
#            #genome=os.getenv('1')#link to $1 in batch script
#            #genome=str(file)+".fa"
#            genome='fulGla.fa'
#            z = SeqIO.index_db(genome+".idx", genome, "fasta")
#            z = z.get(tName)[tStart:tEnd] #returns the target exon from the target scaffold as a SeqRecord 
#            fasta = z.format("fasta") #returns the SeqRecord in fasta format
#            if rc:
#                zrc=(z.reverse_complement(id="rc_",description=z.description))
#                rcfasta=zrc.format("fasta")
#                fout.write(rcfasta)
#            else:
#                fout.write(fasta)
#    fout.close()
