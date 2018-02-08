# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 16:50:00 2015

@author: Phil
"""

"""
Created on Wed Sep 16 15:28:34 2015

@author: phil-grayson
"""
import os 
import csv
import re
if not os.path.exists('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds'): #if this path doesn't exist
    os.makedirs('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds') #make it
           
i = 1 #i is file number
j = 1 #j is folder number
x = [] #empty list for matches
y = [] #empty list to start collecting bed lines in

with open('/n/home12/pgrayson/regal/PseudoSearch/Final/Chicken_CDS.bed', 'rU') as handle: #opens bed in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    if not os.path.exists('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s' % j): #if this path doesn't exist
        os.makedirs('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s' % j) #make it
    for strLine in reader:
        if len(x)==0: #if the program has just started, there is nothing in x
            pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #search for ID=cds#                          
            pcmatch = pcsearch.group(1) #everything in re brackets becomes pcmatch
            x = [pcmatch] #and now x is pcmatch too
            if i < 650: #if the file number in the folder is under 650
                if str(pcmatch) == strLine[3]: #if the current line has pcmatch in column 4
                    y.append(strLine) #append that line to list y 
            else: #once we hit i=650, make a new folder, and place 650 into that folder as i=1
                j += 1 #here we add 1 to j (the folder number)
                os.makedirs('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s' % j)
                i = 1 #reset i (the file number)
                if str(pcmatch) in strLine[3]: #same as above
                    y=[]                    
                    y.append(strLine)
        if x[0] == str(pcmatch):
            if i < 650:
                if str(pcmatch) in strLine[3]: #same as above
                    y.append(strLine)
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
                j += 1
                os.makedirs('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s' % j)
                i = 1
                if str(pcmatch) in strLine[3]: #same as above
                    y.append(strLine)
                if str(pcmatch) not in strLine[3]: #if pcmatch is not in column 4 of a line (i.e., we have all exons to a transcript)
                    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #write a file
                    fout= open(fileout, 'w')
                    while len(y)>0: #since y is a list of lists (one list for each exon in the transcript) we pop it until it is len==0
                        y_item = (y.pop(0)) #pop the first list (bed line) off
                        y_item_tab = '\t'.join(y_item) #tab seperate it
                        fout.write(y_item_tab + '\n') #and write it to the file
                    fout.close() #close that file.  we're done with it
                    y = [] #start y over
                    i += 1 #add 1 to the file counter
                    y.append(strLine) #add the first line of the new transcript to y
                    pcsearch = re.search(r'(ID\=cds\d*)\;', strLine[3]) #same as above           
                    pcmatch = pcsearch.group(1)
                    x = [pcmatch]
    fileout=('/n/home12/pgrayson/regal/PseudoSearch/Final/Beds/%s/%s_%s.bed' % (j,j,i)) #here because after every line has been read, we haven't popped y into a file
    fout= open(fileout, 'w')
    while len(y)>0: 
        y_item = (y.pop(0)) 
        y_item_tab = '\t'.join(y_item)
        fout.write(y_item_tab + '\n') 
    fout.close()
