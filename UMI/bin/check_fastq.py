#校验fastq文件
from  sys import argv
import gzip
import re

infile = argv[1]

def check_orchard(seq):
    for s in seq.strip("\n"):
        if s not in ["A","C","G","T","N"]:
            return False
    return True

t = 0
for line in gzip.open(infile):
    if t % 4 == 0:
        if not re.search("\/1|\/2", str(line, encoding="utf-8").strip("\n")):
            print(line)
    if t % 4 == 1:
        bool = check_orchard(str(line, encoding="utf-8"))
        if not bool:
            print(line)
    t += 1
