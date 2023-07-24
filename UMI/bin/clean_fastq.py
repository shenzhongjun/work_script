"""
处理fastq文件，如果第7位碱基不是T，则去除
"""
import gzip
import re
from argparse import ArgumentParser
import os
from multiprocessing import Pool
from Bio import SeqIO

class clean_data:
    def __init__(self, infile, outfile, cpu, tempdir):
        self.infile = infile
        self.outfile = outfile
        self.cpu = cpu
        self.tempdir = tempdir

    def judge(self, record):
        seq = record.seq
        if seq[6] != "T":
            new_seq = "NNNNNNN" + seq[7:]
        else:
            new_seq = seq
        return str(new_seq)

    def quality_control(self, infile, outfile):
        f = gzip.open(infile, "rt")
        o = open(outfile, "w")
        for record in SeqIO.parse(f, 'fastq'):
            read_name = str(record.id)
            new_seq = self.judge(record)
            quality = record.letter_annotations['phred_quality']
            ascii_qulity = [chr(int(i) + 33) for i in quality]
            #if len(ascii_qulity) > len(new_seq):
            #    ascii_qulity = ascii_qulity[:len(new_seq)]
            outline = ["@" + read_name, new_seq, "+", "".join(ascii_qulity)]
            o.write("\n".join(outline)+"\n")
        o.close()
        os.system(f"rm {infile}")

    def compress(self, file):
        cmds = f"gzip -c {file} > {file}.gz"
        os.system(cmds)
        cmds = f"rm {file}"
        os.system(cmds)

    def deploy(self):
        os.makedirs(self.tempdir, exist_ok=True)
        cmds = f"zcat {self.infile} | split -l 20000000 -d  -  {self.tempdir}/part_ "
        os.system(cmds)
        cmds = fr"for i in $(ls {self.tempdir}/part_*);do mv $i $i.fastq;done"
        os.system(cmds)

        pool = Pool(self.cpu)
        for f in os.listdir(self.tempdir):
            if not re.search("fastq$", f):
                continue
            pool.apply_async(self.compress, (os.path.join(self.tempdir, f), ))
        pool.close()
        pool.join()

        pool = Pool(self.cpu)
        result = []
        for f in os.listdir(self.tempdir):
            if not re.search("fastq.gz$", f):
                continue
            infile = os.path.join(self.tempdir, f)
            outfile = os.path.join(self.tempdir, f.replace(".fastq.gz", ".clean.fastq.gz"))
            result.append(pool.apply_async(self.quality_control, (infile, outfile, )))
        pool.close()
        pool.join()

        files = os.listdir(self.tempdir)
        all_file = " ".join([os.path.join(self.tempdir, f) for f in files])
        outfile = os.path.join(self.tempdir, "merge.fastq")
        cmds = f"cat {all_file} >> {outfile}"
        os.system(cmds)
        cmds = f"gzip -c {outfile} > {self.outfile}"
        os.system(cmds)
        cmds = f"rm {self.tempdir}/*"
        os.system(cmds)
def main():
    parser = ArgumentParser()
    parser.add_argument("-i", action="store", required=True, help="input fastq file")
    parser.add_argument("-o", action="store", required=True, help="output fastq file")
    parser.add_argument("--cpu", action="store", default=8, type=int, help="cpu")
    parser.add_argument("--tempdir", action="store", required=True, help="tempdir")
    argv = parser.parse_args()
    t = clean_data(argv.i, argv.o, argv.cpu, argv.tempdir)
    t.deploy()

if __name__ == '__main__':
    main()