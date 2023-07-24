# /usr/bin/env python
# conding:utf-8
#从bam中挑选reads
import argparse
import pysam
import subprocess
import random
import os
from multiprocessing import Pool

def countRead(bam, samtools):
    status, output = subprocess.getstatusoutput("{samtools} view -c {bam}".format(bam=bam, samtools=samtools))
    return output

class pickBam:
    def __init__(self, total_Duplex, bam, outDir, cpu, samtools):
        self.total_Duplex = total_Duplex
        self.bam = bam
        self.outDir = outDir
        self.cpu = cpu
        self.samtools = samtools

    def pickNum(self):
        numbers = []
        for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            numbers.append(int(int(total_Duplex) * i))
        self.numbers = numbers
        print(self.numbers)

    def write_read(self, order, outfile):
        bamFile = pysam.AlignmentFile(self.bam, "rb")
        out = pysam.AlignmentFile(outfile, "wb", template=bamFile)
        j = 0
        for read in bamFile.fetch():
            r = random.random()
            if r > (order + 1) * 0.1 + 0.05:
                continue
            if j > self.numbers[order]:
                print(str(order) + "\t" + str(j))
                break
            j += 1
            out.write(read)
        out.close()

    def pickBam(self):
        outfile = []
        for i in range(len(self.numbers)):
            outfile.append(os.path.join(self.outDir, "pick_{n}.bam".format(n=i)))

        pool = Pool(self.cpu)
        for i in range(len(self.numbers)):
            pool.apply_async(self.write_read, (i, outfile[i], ))
        pool.close()
        pool.join()

    def run(self):
        self.pickNum()
        self.pickBam()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="UMI qc")
    parser.add_argument('--input', help="input bam file", required=True)
    parser.add_argument('--outDir', help="output bam directory file", required=True)
    parser.add_argument('--samtools', help="samtools", default="/mnt/share01/tools/bin/samtools")
    parser.add_argument('--cpu', help="cpu number", type=int, default=8)

    argv = parser.parse_args()
    total_Duplex = countRead(argv.input, argv.samtools)
    print(total_Duplex)
    t = pickBam(total_Duplex, argv.input, argv.outDir, argv.cpu, argv.samtools)
    t.run()