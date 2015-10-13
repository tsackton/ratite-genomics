"""
created oct7-15
@author: phil-grayson
SeqIO.index_db creates a fasta index file.  This was passed over each genome.
"""


from Bio import SeqIO
import os
import sys

file = os.getenv('FILE')
z = SeqIO.index_db(file+".idx", file, "fasta")
