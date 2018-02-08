import sys
import re
import csv
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
exoncount=0
phile = str(sys.argv[1])

class Vividict(dict): #used for autovivification of main dictionary 
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

class Exon(dict): #this class is utilized with the galGal sequences
    def __init__(self, refStart, refEnd): 
        self.refStart = refStart
        self.refEnd = refEnd
        self.qStart = None
        self.qEnd = None
    
    def reflen(self): #this class is also used with galGal
        exonlen=(self.refEnd - self.refStart)        
        return exonlen
    
    def leftPadStr(self): #this class compares exon.refStart and exon.refEnd to the qStart and qEnd of each individual exon in the target genome to produce leftPads
        nlen=abs(int(self.refStart - self.qStart))
        return("N"*nlen)
         
    def rightPadStr(self): #as above, but for rightPads 
        nlen=abs(int(self.refEnd - self.qEnd))
        return("N"*nlen)     


try:
    outlist = []    
    with open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/dic/'+phile.split('.')[0]+".txt",'rU') as outlog:
        for eachline in outlog:
            outlist.append(eachline.strip('\n'))
except:
    outlist = []
            
        
cdic = Vividict() #the main dictionary
exonkey={} #the dictionary to keep track of exon number
with open('/n/home12/pgrayson/regal/PseudoSearch/Final/HLO_psl/singles/singles_galGal.psl', 'rU') as handle: #opens chicken psl (our reference) in universal mode
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #line by line, calls the first column name
        name = strLine[0]        
        try:
            namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',name) #grab ###,Genbank:??_####
            rStart = int(strLine[16])
            rEnd = int(strLine[17])
            if namegrab.group(1) not in cdic:
                exoncount = 1 #begins exon count when the namegrab group is not already in cdic
            else:
                exoncount +=1 #if namegrab group is in cdic, add to exon count
            cdic[namegrab.group(1)][int(exoncount)] = Exon(rStart,rEnd) #this line places all important info into cdic.  the namegrab group, the exon #, and ref start and end from the Exon class
            exonkey[name]=int(exoncount) #this line gives us our exon count dictionary with the entire name field as key and exon # as the value
        except:
            None #necessary to eliminate CDS entries that do not code for proteins


with open('/n/home12/pgrayson/regal/PseudoSearch/Final/HLO_psl/singles/singles_'+phile.split('.')[0]+'.psl', 'rU') as handle: #opens the target genome
    reader=csv.reader(handle,delimiter='\t') #reads file with tabs as delimiters
    for strLine in reader: #read the tab delimited file in and identify the necessary columns
        name = strLine[0]
        tStart = int(strLine[16])
        tEnd = int(strLine[17]) 
        tName = strLine[14]
        qStart = int(strLine[12])
        qEnd = int(strLine[13])
        strand = strLine[9]        
        try: 
            namegrab = re.search('[a-zA-z\/\-]*\:([P0-9]*\,Genbank\:[A-Za-z0-9/./_]*)\,(.*\:[\+\-])',name) #as above
            exnum = int(exonkey[name]) #for each exon, we pair it to its correct exon number from chicken using our exonkey dict
            cdic[namegrab.group(1)][exnum].update({"tStart":tStart,"tEnd":tEnd,"Scaff":tName,"Strand":strand}) #add all of the target info to cdic.  Note that this call differs from above in that it does not utilize the Exon class.  This allows us to identify missing exons in our target further down the script.  It also includes much more information (strand and scaffold)
            cdic[namegrab.group(1)][exnum].qStart=qStart #assign chicken qStart for given exon for use in Exon class to determine pads
            cdic[namegrab.group(1)][exnum].qEnd=qEnd #assign chicken qEnd for given exon as above
        except:
            None #this is necessary to allow the program to continue past exons that are missing in our target


try:
    flog = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/dic/'+phile.split('.')[0]+".txt",'a')
    fout = open('/n/home12/pgrayson/regal/PseudoSearch/Final/HLO_psl/singles/outSingles/concat_'+phile,'a')
except:
    flog = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/dic/'+phile.split('.')[0]+".txt",'w')
    fout = open('/n/home12/pgrayson/regal/PseudoSearch/Final/HLO_psl/singles/outSingles/concat_'+phile,'w') #multifasta out file is created
exonlist = [] #exonlist is used to collect all exons from a given transcript
genome=phile #call the genome
gindex = SeqIO.index_db(genome+".idx", genome, "fasta") #call to the genome and its index
for trans in cdic.keys(): #for a given transcript
    if str(trans) not in outlist:
        st = 'None'
        description = 'None'
        for exon in sorted(cdic[trans].keys()): #for each exon of that transcript
            if len (cdic[trans][exon].keys()) == 0: #if the exon length is 0
                missingexon = str("N"*cdic[trans][exon].reflen()) #missingexon is created as a string of N's that equals the exon length in chicken
                z = SeqRecord(Seq(missingexon)) #create a SeqRecord for this string of N's
                exonlist.append(z.seq) #append this SR to exonlist
            else: #if the exon is not 0-length
                z = gindex.get(cdic[trans][exon]['Scaff'])[cdic[trans][exon]['tStart']:cdic[trans][exon]['tEnd']] #pull the full length of aligned target exon from the target genome
                z.seq=(cdic[trans][exon].leftPadStr())+z.seq+(cdic[trans][exon].rightPadStr()) #place the appropriate # of N's onto each end of this sequence to match it to the exon length in chicken
                if cdic[trans][exon]['Strand'] == "+-": #here is where we introduce reverse complement as a variable based on the psl's strand column in our cdic
                    rc = True
                elif cdic[trans][exon]['Strand'] == "--":
                    rc = True
                else:
                    rc = False
                description=z.description #save the description of the SeqRecord that has been pulled (used in final fasta header)
                st=str(cdic[trans][exon]['Strand']) #create st variable to keep track of strand in final fasta header
                if rc: #if rc is True, reverse complement the sequence
                    zrc=(z.reverse_complement(description=description))
                    exonlist.append(zrc.seq) #and append it to the exonlist
                else: #if not, just append the exon
                    exonlist.append(z.seq)
        transcript="".join([str(seq_rec) for seq_rec in exonlist]) #join all the exons in exonlist in the right order    
        exonlist=[] #empty the list
        fout.write(">ChickenGeneID:"+trans+",pslStrand:"+st+","+description+"\n"+transcript+"\n") #write them out with appropriate label
        fout.flush() #flush internal buffer
        os.fsync(fout.fileno()) #flush os buffer
        flog.write(trans+"\n")
        flog.flush()
        os.fsync(flog.fileno())
       
fout.close() #close the outfile
flog.close() #close the logfile
