#校验fastq文件
from sys import argv
import gzip
import re
from Bio import SeqIO


infile = argv[1]
read = argv[2]
f = gzip.open(infile, "rt")
for record in SeqIO.parse(f, 'fastq'):
    read_name = str(record.id)
    if read in read_name:
        print(record.seq)





