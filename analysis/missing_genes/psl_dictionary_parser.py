import sys
import os
import re
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
exoncount=0
os.chdir('/Users/Phil/Desktop/')

class Vividict(dict): #used as autovivification
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class Exon(dict): #this class 
    def __init__(self, refStart, refEnd): 
        self.refStart = refStart
        self.refEnd = refEnd
        self.qStart = None
        self.qEnd = None
    
    def reflen(self):
        exonlen=(self.refEnd - self.refStart)        
        return exonlen
    
    def leftPadStr(self):
        nlen=abs(int(self.refStart - self.qStart))
        return("N"*nlen)
         
    def rightPadStr(self):
        nlen=abs(int(self.refEnd - self.qEnd))
        return("N"*nlen)     
        
cdic = Vividict()
exonkey={}
with open('/Users/Phil/Desktop/singles_galGal.psl', 'rU') as handle: #opens psl in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in call the first column name
        name = strLine[0]        
        try: 
            namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',name)
            rStart = int(strLine[16])
            rEnd = int(strLine[17])
            if namegrab.group(1) not in cdic:
                exoncount = 1
            else:
                exoncount +=1
            cdic[namegrab.group(1)][int(exoncount)] = Exon(rStart,rEnd)
            exonkey[name]=int(exoncount)
        except:
            None #necessary for CDS entries that do not code for proteins


with open('/Users/Phil/Desktop/fulGla_singles.txt', 'rU') as handle:
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in identify the necessary columns
        name = strLine[0]
        tStart = int(strLine[16])
        tEnd = int(strLine[17]) 
        tName = strLine[14]
        qStart = int(strLine[12])
        qEnd = int(strLine[13])
        strand = strLine[9]        
        try: 
            namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',name)
            exnum = int(exonkey[name])
            cdic[namegrab.group(1)][exnum].update({"tStart":tStart,"tEnd":tEnd,"Scaff":tName,"Strand":strand}) #add all of the target info to cdic
            cdic[namegrab.group(1)][exnum].qStart=qStart #assign chicken qStart for that exon
            cdic[namegrab.group(1)][exnum].qEnd=qEnd #assign chicken qEnd for that exon
        except:
            None

fout = open('/Users/Phil/Desktop/fulGla_concat.txt','w') #multifasta file
exonlist = []
genome='fulGla.fa'      
for trans in cdic.keys():
    for exon in sorted(cdic[trans].keys()):
        if len (cdic[trans][exon].keys()) == 0:
            missingexon = str("N"*cdic[trans][exon].reflen())
            z = SeqRecord(Seq(missingexon))
            exonlist.append(z.seq)
        else:
            z = SeqIO.index_db(genome+".idx", genome, "fasta")
            z = z.get(cdic[trans][exon]['Scaff'])[cdic[trans][exon]['tStart']:cdic[trans][exon]['tEnd']]
            z.seq=(cdic[trans][exon].leftPadStr())+z.seq+(cdic[trans][exon].rightPadStr())
            #exonlist.append(z.seq)
            if cdic[trans][exon]['Strand'] == "+-":
                rc = True
            elif cdic[trans][exon]['Strand'] == "--":
                rc = True
            else:
                rc = False
            description=z.description
            if rc:
                zrc=(z.reverse_complement(description=description))
                exonlist.append(zrc.seq)
            else:
                exonlist.append(z.seq)
    transcript="".join([str(seq_rec) for seq_rec in exonlist])    
    exonlist=[]
    if rc:
        fout.write(">"+"RC_"+description+"\n"+transcript+"\n") #write them out with appropriate label
    else:
        fout.write(">"+description+"\n"+transcript+"\n") #as above   
fout.close()