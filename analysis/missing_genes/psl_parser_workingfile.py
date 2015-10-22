import csv
from Bio import SeqIO
import os
import re
x = [] #empty list for gene ID matches
y = [] #empty list to start collecting psl lines in  

#fout = open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/counts/'+os.getenv('FILE')+'.text','w') #out files will be ${FILE}.txt
#with open('/n/regal/edwards_lab/phil/PseudoSearch/Final/HLO_psl/'+os.getenv('FILE'), 'rU') as handle:
fout = open('/Users/Phil/Desktop/messingaroundfulGla.txt','w') #multifasta file
fout2 = open('/Users/Phil/Desktop/messingaroundfulGla2.txt','w') #table with exon numbers and classifications of exons
os.chdir('/Users/Phil/Desktop/')
with open('/Users/Phil/Desktop/1fulGla.psl', 'rU') as handle: #opens psl in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in call the first column name
        name = strLine[0]
        if len(x)==0:  #if the program has just started, there is nothing in x      
            namegrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,.*\:[\+\-]',name) #grab (gene ID) from name                    
            genegrab = namegrab.group(1) #group 1 becomes genegrab          
            x = [genegrab] #and now x is genegrab too
            if genegrab in name: #if genegrab is in name
                y.append(strLine) #append the whole line to y
        else: #if the script hasn't just started
            nextgrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,.*\:[\+\-]',name) #as above, but called nextgrab 
            if x[0] == nextgrab.group(1): #if x[0] from above matches nextgrab              
                y.append(strLine) #append the line
            else: #if they don't match, we have all exons for a transcript             
                total = 0 #a bunch of counts that start at 0
                perfect = 0
                imperfect = 0
                falseStart = 0
                falseEnd = 0
                falseBoth = 0                            
                for listLine in y: #y is a list of psl lines for a transcript, label the columns
                    name = listLine[0]                    
                    strand = listLine[9]
                    qName = listLine[10]
                    qSize = int(listLine[11])
                    qStart = int(listLine[12])
                    qEnd = int(listLine[13])
                    tName = listLine[14]
                    tSize = int(listLine[15])
                    tStart = int(listLine[16])
                    tEnd = int(listLine[17])
                    blockCount = int(listLine[18])
                    blockSizes = listLine[19]
                    qStarts = listLine[20]
                    tStarts = listLine[21]
                    qnamestart = re.search('Start\:([0-9]*)\,', name) #pull the chicken start position from name
                    qnamestop = re.search('Stop\:([0-9]*)\,', name) #pull the chicken end position from name
                    if int(qnamestart.group(1)) == qStart and int(qnamestop.group(1)) == qEnd: #if the chicken positions in name match those used in qStart and qEnd
                        qExonSize=qEnd-qStart #check qExon size
                        tExonSize=tEnd-tStart #check tExon size
                        if qExonSize==tExonSize: #if both these exons are the same size
                            perfect +=1 #this is a perfect sequence
                            listLine.append("perfect") #add "perfect" to the listLine for later
                            listLine.append("0")
                            listLine.append("0")
                        else:
                            imperfect +=1 #otherwise, it is imperfect (i.e., it has indels)
                            listLine.append("imperfect")
                            listLine.append("0")
                            listLine.append("0")
                    else: #if the entire chicken exon is not captured, there are 3 more categories
                        if int(qnamestart.group(1)) != qStart and int(qnamestop.group(1)) == qEnd: #missing start, but end is correct
                            falseStart +=1
                            fS = abs(int(qnamestart.group(1)) - qStart)
                            listLine.append("falseStart")
                            listLine.append(str(fS))
                            listLine.append("0")
                        elif int(qnamestop.group(1)) != qEnd and int(qnamestart.group(1)) == qStart: #missing end, but start is correct
                            falseEnd +=1
                            fE = abs(int(qnamestop.group(1)) - qEnd)                            
                            listLine.append("falseEnd")
                            listLine.append("0")
                            listLine.append(str(fE))
                        elif int(qnamestop.group(1)) != qEnd and int(qnamestart.group(1)) != qStart: #missing start and end
                            falseBoth +=1
                            fS = abs(int(qnamestart.group(1)) - qStart)
                            fE = abs(int(qnamestop.group(1)) - qEnd)
                            listLine.append("falseBoth")
                            listLine.append(str(fS))
                            listLine.append(str(fE))
                if falseStart == falseEnd == falseBoth == 0 and perfect + imperfect >= 1: #if y contains a gene with only perfect and imperfect exons                   
                    exonlist = [] #create a list for exons                   
                    exoncount = 0 #start counting exons at 0                   
                    for listLine in y:
                        name = listLine[0] 
                        tName = listLine[14]
                        strand = listLine[9]
                        tStart = int(listLine[16])
                        tEnd = int(listLine[17])
                        classif = str(listLine[22])
                        if strand == "+-": #make sure that reverse compliment is considered
                            rc = True
                        elif strand == "--":
                            rc = True
                        else:
                            rc = False
                        genome='fulGla.fa' #import the genome
                        z = SeqIO.index_db(genome+".idx", genome, "fasta") #call the genome's index
                        z = z.get(tName)[tStart:tEnd] #return the target exon from the target scaffold as a SeqRecord 
                        if rc:
                            zrc=(z.reverse_complement(id="rc_",description=z.description)) #rc the piece
                            exonlist.append(zrc.seq) #add it to exon list
                            exoncount +=1 #add to exon count
                            listLine.append("Exon_"+str(exoncount)) #append the exon number to listLine
                            listLine2 = '\t'.join(listLine) #turn that line (which is in list format) into a tab delimited string                        
                            fout2.write(listLine2 + '\n') #write that line to the table file
                        else:
                            exonlist.append(z.seq) #just normal, so no rc
                            exoncount +=1 #as above
                            listLine.append("Exon_"+str(exoncount)) 
                            listLine2 = '\t'.join(listLine)                        
                            fout2.write(listLine2 + '\n')
                    transcript="".join([str(seq_rec) for seq_rec in exonlist]) #concatenate all exons in exonlist
                    if rc:
                        fout.write(">"+"RC_"+z.description+"\n"+transcript+"\n") #write them out with appropriate label
                        y=[]
                        y.append(strLine)
                    else:
                        fout.write(">"+z.description+"\n"+transcript+"\n") #as above
                        y=[]
                        y.append(strLine)
                else: #if we do not have only perfect and imperfect exons in the transcript                   
                    exonlist = [] #all as above                   
                    exoncount = 0                   
                    for listLine in y:
                        name = listLine[0] 
                        tName = listLine[14]
                        strand = listLine[9]
                        tStart = int(listLine[16])
                        tEnd = int(listLine[17])
                        classif = str(listLine[22])
                        falS = int(listLine[23])
                        falE = int(listLine[24])
                        if strand == "+-":
                            rc = True
                        elif strand == "--":
                            rc = True
                        else:
                            rc = False
                        genome='fulGla.fa'
                        z = SeqIO.index_db(genome+".idx", genome, "fasta")
                        z = z.get(tName)[tStart:tEnd]
                        if "perfect" in classif: #if a given exon is perfect or imperfect (all as above)                            
                            if rc: 
                                zrc=(z.reverse_complement(id="rc_",description=z.description))
                                exonlist.append(zrc.seq)
                                exoncount +=1
                                listLine.append("Exon_"+str(exoncount))
                                listLine2 = '\t'.join(listLine)                        
                                fout2.write(listLine2 + '\n')
                            else:
                                exonlist.append(z.seq)
                                exoncount +=1
                                listLine.append("Exon_"+str(exoncount)) 
                                listLine2 = '\t'.join(listLine)                    
                                fout2.write(listLine2 + '\n')
                        if "perfect" not in classif:
                            if classif == "falseStart": #if the exon isn't perfect or imperfect, we classify and add ns as appropriate
                                z.seq=("N"*falS)+z.seq
                            if classif == "falseEnd":
                                z.seq=z.seq+("N"*falE)
                            if classif == "falseBoth":
                                z.seq=("N"*falS)+z.seq+("N"*falE)
                            if rc: 
                                zrc=(z.reverse_complement(id="rc_",description=z.description))
                                exonlist.append(zrc.seq)
                                exoncount +=1
                                listLine.append("Exon_"+str(exoncount))
                                listLine2 = '\t'.join(listLine)                        
                                fout2.write(listLine2 + '\n')
                            else:
                                exonlist.append(z.seq)
                                exoncount +=1
                                listLine.append("Exon_"+str(exoncount)) 
                                listLine2 = '\t'.join(listLine)                    
                                fout2.write(listLine2 + '\n')            
                    transcript="".join([str(seq_rec) for seq_rec in exonlist]) #concatenate all exons in exonlist
                    if rc:
                        fout.write(">"+"RC_"+z.description+"\n"+transcript+"\n") #write them out with appropriate label
                        y=[]
                        y.append(strLine)
                    else:
                        fout.write(">"+z.description+"\n"+transcript+"\n") #as above
                        y=[]
                        y.append(strLine)
            x =[nextgrab.group(1)]
fout.close()
fout2.close()                    

#                            
                
                
    
        
        
        
        if str(pcmatch) not in strLine[3]: #if pcmatch is not in column 4 of a line (i.e., we have all exons to a transcript)
                    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #write a file
                    fout= open(fileout, 'w')
                    while len(y)>0: #since y is a list of lists (one list for each exon in the transcript) we pop it until it is len==0
                        y_item = (y.pop(0)) #pop the first list (bed line) off
                        y_item_tab = '\t'.join(y_item) #tab seperate it
                        fout.write(y_item_tab + '\n') #and write it to the file
                    fout.close() #close that file -  we're done with it
                    y = [] #start y over
                    i += 1 #add 1 to the file counter
                    y.append(strLine) #add the first line of the new transcript to y
                    pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #same as above           
                    pcmatch = pcsearch.group(1)
                    x = [pcmatch]    
    fout.close()
           
                if str(pcmatch) not in strLine[3]: #if pcmatch is not in column 4 of a line (i.e., we have all exons to a transcript)
                    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #write a file
                    fout= open(fileout, 'w')
                    while len(y)>0: #since y is a list of lists (one list for each exon in the transcript) we pop it until it is len==0
                        y_item = (y.pop(0)) #pop the first list (bed line) off
                        y_item_tab = '\t'.join(y_item) #tab seperate it
                        fout.write(y_item_tab + '\n') #and write it to the file
                    fout.close() #close that file -  we're done with it
                    y = [] #start y over
                    i += 1 #add 1 to the file counter
                    y.append(strLine) #add the first line of the new transcript to y
                    pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #same as above           
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

            
            qExonSize=qEnd-qStart
            tExonSize=tEnd-tStart
            if qExonSize==tExonSize:
                #genome=os.getenv('1')#link to $1 in batch script
                #genome=str(file)+".fa"
                genome='fulGla.fa'
                z = SeqIO.index_db(genome+".idx", genome, "fasta")
                z = z.get(tName)[tStart:tEnd] #returns the target exon from the target scaffold as a SeqRecord 
                fasta = z.format("fasta") #returns the SeqRecord in fasta format
                if rc:
                    zrc=(z.reverse_complement(id="rc_",description=z.description))
                    rcfasta=zrc.format("fasta")
                    fout.write(rcfasta)
                else:
                    fout.write(fasta)
    fout.close()
            
        
        else:
            #here is where we write the n's in.
        if int(qnamestart.group(1)) != qStart:
            print "Start is off for "+str(strLine)
        if int(qnamestop.group(1)) != qEnd:
            print "End is off for "+str(strLine)
        
        if strand == "+-":
            rc = True
        elif strand == "--":
            rc = True
        else:
            rc = False
        qExonSize=qEnd-qStart
        tExonSize=tEnd-tStart
        if qExonSize==tExonSize:
            #genome=os.getenv('1')#link to $1 in batch script
            #genome=str(file)+".fa"
            genome='fulGla.fa'
            z = SeqIO.index_db(genome+".idx", genome, "fasta")
            z = z.get(tName)[tStart:tEnd] #returns the target exon from the target scaffold as a SeqRecord 
            fasta = z.format("fasta") #returns the SeqRecord in fasta format
            if rc:
                zrc=(z.reverse_complement(id="rc_",description=z.description))
                rcfasta=zrc.format("fasta")
                fout.write(rcfasta)
            else:
                fout.write(fasta)
    fout.close()


##if we need to do different things for different rev comps :
#                            exontot = len(y)+1
#                            exoncount += 1
#                            listLine.append("Exon_"+str(exontot-exoncount))