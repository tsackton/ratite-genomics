#will grab the translated concatenated mRNA transcript from every .fa file in a given directory.  Input is just chickenID number

from Bio import SeqIO
import re
import sys
import os

fout = open('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/mfr/'+str(sys.argv[1])+".fa",'w')
for filename in os.listdir('/n/home12/pgrayson/regal/PseudoSearch/genomes/raw/concat/translation'):
    if filename.endswith('_transbestframe.txt'):
        for record in SeqIO.parse(filename,"fasta"):
            recordgrab = re.search('ChickenGeneID\:([0-9]*.*Genbank\:[A-Za-z0-9/./_]*)',record.description)
            if recordgrab.group(1).split(",")[0] == str(sys.argv[1]).split(",")[0]:
                trans = record.seq
                fout.write(">"+filename.split("_")[0]+"_"+str(record.description)+"\n"+str(trans)+"\n")
fout.close()
