# /usr/bin/env python
# conding:utf-8
#插入片段长度质控
from argparse import ArgumentParser
import pysam
from multiprocessing import Pool
import matplotlib.pyplot as plt

class count_insert_size:
    def __init__(self, bam, bed, cpu, outfile):
        self.bam = bam
        self.bed = bed
        self.cpu = cpu
        self.outfile = outfile

    def count_depth(self, line):
        bamFile = pysam.AlignmentFile(self.bam, "rb")
        InsertSizeDict = {}
        lines = line.split("\t")
        ###获取插入片段长度
        for alignment in bamFile.fetch(lines[0], int(lines[1]), int(lines[2]), until_eof=True):
            if alignment.is_unmapped:
                continue
            if alignment.template_length == 0:
                continue
            InsertSizeDict.setdefault(alignment.query_name, int(abs(alignment.template_length)))
        bamFile.close()
        return InsertSizeDict

    def merge_dict(self, InsertSizeDict, Length_Dict):
        for readID, length in InsertSizeDict.items():
            if length in Length_Dict.keys():
                Length_Dict[length] += 1
            else:
                Length_Dict[length] = 1
        return Length_Dict

    def write_file(self, Length_Dict):
        out = open(self.outfile, "w")
        length_key = sorted(Length_Dict.keys())
        for line in length_key:
            out.write(str(line) + "\t" + str(Length_Dict[line]) + "\n")
        out.close()

    def deploy(self):
        pool = Pool(self.cpu)
        result = []
        for line in open(self.bed, "r").readlines():
            line = line.strip()
            result.append(pool.apply_async(self.count_depth, (line, )))
        pool.close()
        pool.join()
        Length_Dict = {}
        for res in result:
            InsertSizeDict = res.get()
            Length_Dict = self.merge_dict(InsertSizeDict, Length_Dict)
        self.write_file(Length_Dict)

class inster_size_qc:
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile
        self.cutoff = 0.8

    def qc(self, Length_dict):
        total_size = 0
        size_below_250 = 0
        length_key = sorted(Length_dict.keys())
        for line in length_key:
            total_size += 1
            if line != 0 and Length_dict[line] <= 250:
                size_below_250 += 1
        size_below_250_ratio = round(size_below_250 / total_size, 4)
        return size_below_250_ratio

    def judge_status(self, ratio):
        status = ""
        if ratio < self.cutoff:
            status = "Fail"
        else:
            status = "Pass"
        return status

    def write(self, size_below_250_ratio, status, Peak_Instert_Size, outfile):
        out = open(outfile, "w")
        out.write("Peak_Instert_Size\tPCT(250bp)\tStatus\n")
        out.write(str(Peak_Instert_Size)+"\t"+str(size_below_250_ratio) + "\t" + status+"\n")
        out.close()

    def peak_Size(self, Length_dict):
        Peak_Instert_Size = 0
        Max_Freq_Num = 0
        length_key = sorted(Length_dict.keys())
        for line in length_key:
            if int(line) < 10:
                continue
            if int(Length_dict[line]) > Max_Freq_Num:
                Peak_Instert_Size = int(line)
                Max_Freq_Num = int(Length_dict[line])
            else:
                pass
        return Peak_Instert_Size

    def deploy(self):
        Length_dict = {}
        for line in open(self.infile, "r").readlines():
            Length, num = line.strip().split("\t")
            Length_dict.setdefault(int(Length), int(num))
        size_below_250_ratio = self.qc(Length_dict)
        status = self.judge_status(size_below_250_ratio)
        Peak_Instert_Size = self.peak_Size(Length_dict)
        self.write(size_below_250_ratio, status, Peak_Instert_Size, self.outfile)

def get_data(infile):
    xdata = []
    ydata = []
    whole = 0
    for line in open(infile, "r").readlines():
        lineinfo = line.strip('\n').split('\t')
        whole += int(lineinfo[1])
        if int(lineinfo[0]) > 500:
            continue
        xdata.append(int(lineinfo[0]))
        ydata.append(int(lineinfo[1]))
    return xdata, ydata

def plot_curve(file1, file2, file3, outfile, sample, label1, label2):
    #print(whole)
    #print(ydata[0])
    xdata1, ydata1 = get_data(file1)
    xdata2, ydata2 = get_data(file2)
    xdata3, ydata3 = get_data(file3)
    plt.figure(figsize=(20, 8))
    #plt.bar(range(len(xdata)), ydata)
    plt.plot(range(len(xdata1)), ydata1, label=label1)
    plt.plot(range(len(xdata2)), ydata2, label=label2)
    plt.plot(range(len(xdata3)), ydata3, label="3")
    plt.xlabel("InsertSize")
    plt.legend(fontsize=20)
    plt.ylabel("Count")
    plt.title(sample + " Insert Size Distribution Plot")
    plt.savefig(outfile)
    plt.close()

if __name__ == '__main__':
    parser = ArgumentParser()
    #parser.add_argument("--bam", required=True, help="bamfile after markdup")
    #parser.add_argument("--bed", required=True, help="bed file")
    parser.add_argument("--sample", required=True, help="sample name")
    parser.add_argument("--file1", required=True, help="output InsertSize file")
    parser.add_argument("--file2", required=True, help="output InsertSize file")
    parser.add_argument("--file3", required=True, help="output InsertSize file")
    #parser.add_argument("--outfile", required=True, help="output qc status file")
    parser.add_argument("--png", required=True, help="output png file")
    parser.add_argument("--label1", required=True)
    parser.add_argument("--label2", required=True)
    #parser.add_argument("--cpu", action="store", default=8, help="number of cpus")
    o = parser.parse_args()
    #统计插入片段，并将结果写入文档
    #t = count_insert_size(o.bam, o.bed, o.cpu, o.InsertSize)
    #t.deploy()
    #对插入片段的分布进行质控
    #t = inster_size_qc(o.InsertSize, o.outfile)
    #t.deploy()
    #画图
    plot_curve(o.file1, o.file2, o.file3, o.png, o.sample, o.label1, o.label2)
