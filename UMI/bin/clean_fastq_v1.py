"""
处理fastq文件，如果第7位碱基不是T，则去除
"""
import gzip
from argparse import ArgumentParser
from multiprocessing import Pool
from Bio import SeqIO

class clean_data:
    def __init__(self, infile, outfile, cpu):
        self.infile = infile
        self.outfile = outfile
        self.cpu = cpu

    def run(self, seq, record):
        if seq[6] != "T":
            return False
        else:
            return record

    def deploy(self):
        """
        pool = Pool(self.cpu)
        result = []
        f = gzip.open(self.infile, "rt")
        for record in SeqIO.parse(f, 'fastq'):
            result.append(pool.apply_async(self.run, (record.seq, record, )))
        pool.close()
        pool.join()
        o = gzip.open(self.outfile, "wt")
        for res in result:
            if not res.get():
                continue
            else:
                SeqIO.write(res.get(), o,  "fastq")
        o.close()
        """
        f = gzip.open(self.infile, "rt")
        o = gzip.open(self.outfile, "wt")
        for record in SeqIO.parse(f, 'fastq'):
            result = self.run(record.seq, record)
            if not result:
                continue
            else:
                SeqIO.write(record, o, "fastq")
        o.close()
def main():
    parser = ArgumentParser()
    parser.add_argument("-i", action="store", required=True, help="input fastq file")
    parser.add_argument("-o", action="store", required=True, help="output fastq file")
    parser.add_argument("--cpu", action="store", default=8, type=int, help="cpu")
    argv = parser.parse_args()
    t = clean_data(argv.i, argv.o, argv.cpu)
    t.deploy()

if __name__ == '__main__':
    main()