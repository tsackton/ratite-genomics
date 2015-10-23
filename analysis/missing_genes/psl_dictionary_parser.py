import sys
import os
import re
import csv
exoncount=0
os.chdir('/Users/Phil/Desktop/')

class Vividict(dict): #used as autovivification
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class Exon(dict):
    def __init__(self, Start, End): 
        self.Start = Start
        self.End = End
        self.qStart = None
        self.qEnd = None
    
    def reflen(self):
        exonlen=(self.End - self.Start)        
        return exonlen
    
    def leftPadStr(self):
        nlen=abs(int(self.Start - self.qStart))
        return("N"*nlen)
         
    def rightPadStr(self):
        nlen=abs(int(self.End - self.qEnd))
        return("N"*nlen)     
        
cdic = Vividict()
exonkey={}
with open('/Users/Phil/Desktop/singles_galGal.psl', 'rU') as handle: #opens psl in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in call the first column name
        name = strLine[0]        
        try: 
            namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',name)
            tStart = int(strLine[16])
            tEnd = int(strLine[17])      
            if namegrab.group(1) not in cdic:
                exoncount = 1
            else:
                exoncount +=1
            cdic[namegrab.group(1)][int(exoncount)] = Exon(tStart,tEnd)
            exonkey[name]=int(exoncount)
        except:
            None
        
cdic['425783,Genbank:NP_001264787.1'][1].reflen() #tells how long the exon is

sort(cdic['425783,Genbank:NP_001264787.1'].keys()) #will be useful to sort through exons in order        

"""import sys
import os
import re
import csv
exoncount=0
os.chdir('/Users/Phil/Desktop/')
elist=[]
class Transcript:
    #def __init__(self, name, gene, exon, tStart, tEnd): 
    def __init__(self, name,tStart, tEnd,ename): 
        self.name = name
        self.ename = ename
        #self.gene = gene
        #self.exon = exon
        self.tStart = tStart
        self.tEnd = tEnd
    
    def getName(self):
        return self.name
    
    def __str__(self):
        return "Name: %s, tStart: %s, tEnd: %s, ename: %s" % \
     (self.name, self.tStart, self.tEnd, self.ename)
       # return "Name: %s, Gene: %s, Exon: %s, tStart: %s, tEnd: %s" % \
     #(self.name, self.gene, self.exon, self.tStart, self.tEnd)

ecount=0
with open('/Users/Phil/Desktop/chickentestpsl.txt', 'rU') as handle: #opens psl in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in call the first column name       
        ecount += 1       
        Name = strLine[0]
        namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',Name)
        name = namegrab.group(1)        
        tStart = int(strLine[16])
        tEnd = int(strLine[17])
        #if ename.getName()==namegrab.group(1):
        ename = "line%s" % (ecount)
        ename = Transcript(name,tStart,tEnd,ename)
        elist.append(ename)
    for ename in elist:
        print ename

e.ename


        
print e

fout = open('/Users/Phil/Desktop/messingaroundfulGla.txt','w')
fout.write(str(e.name))
fout.close()



                    strand = listLine[9]
                    qName = listLine[10]
                    qSize = int(listLine[11])
                    qStart = int(listLine[12])
                    qEnd = int(listLine[13])
                    tName = listLine[14]
                    tSize = int(listLine[15])
                    
                    blockCount = int(listLine[18])
                    blockSizes = listLine[19]
                    qStarts = listLine[20]
                    tStarts = listLine[21]
      
        name = strLine[0]
        if len(x)==0:  #if the program has just started, there is nothing in x      
            namegrab = re.search('[a-zA-z\/\-]*\:([0-9]*)\,.*\:[\+\-]',name) #grab (gene ID) from name                    
            genegrab = namegrab.group(1) #group 1 becomes genegrab          
            x = [genegrab] #and now x is genegrab too
            if genegrab in name: #if genegrab is in name
                y.append(strLine) #append the whole line to y"""  