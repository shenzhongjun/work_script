# /usr/bin/env python
# conding:utf-8
#计算一致性双链饱和率
import argparse
import subprocess
import os, re
from multiprocessing import Pool
from statsmodels.nonparametric.smoothers_lowess import lowess
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.interpolate import splev, splrep

def countRead(bam, samtools):
    status, output = subprocess.getstatusoutput("{samtools} view -c {bam}".format(bam=bam, samtools=samtools))
    return output

class lowess_estimate:
    def __init__(self, rbamDir, cbamDir, rbam, cbam, samtools, cpu, png):
        self.rbamDir = rbamDir
        self.cbamDir =cbamDir
        self.rbam = rbam
        self.cbam = cbam
        self.samtools = samtools
        self.cpu = cpu
        self.png = png

    def run_count(self, bam_list):
        count_number = []
        result = []
        pool = Pool(self.cpu)
        for b in bam_list:
            result.append(pool.apply_async(countRead, (b, self.samtools,)))
        pool.close()
        pool.join()
        for res in result:
            count_number.append(res.get())
        count_number.sort()
        return count_number

    def plot(self, xdata, ydata, duplex_1, duplex_2):
        xdata = [int(i) for i in xdata]
        ydata = [int(i) for i in ydata]
        xdata.sort()
        ydata.sort()
        max_reads = max(xdata)
        xdata_d = [round(i / int(max_reads), 2) for i in xdata]
        #plt.scatter(xdata_d, ydata)
        #plt.plot(xdata_d, ydata)
        #平滑曲线
        spl = splrep(xdata_d, ydata)
        x2 = np.linspace(0, 1, 300)
        y2 = splev(x2, spl)

        plt.plot(x2, y2)
        plt.xlim(0, 1)
        plt.ylim(0, duplex_1 * 1.2)
        plt.vlines(0.5, 0, duplex_1, colors="c", linestyles="dashed")
        plt.vlines(1, 0, duplex_1, colors="c", linestyles="dashed")
        plt.hlines(duplex_1, 0, 1, colors="r", linestyles="dashed")
        plt.hlines(duplex_2, 0, 1, colors="r", linestyles="dashed")
        plt.xlabel("Read pairs PCT")
        plt.ylabel("Count of Duplexs")
        plt.title("Duplex yield by read pairs PCT")
        plt.savefig(self.png)
        plt.show()

    def find_bam(self, indir, index):
        bam_list = []
        for f in os.listdir(indir):
            if re.search(index, f):
                bam_list.append(os.path.join(indir, f))
            else:
                pass
        #print(bam_list)
        return bam_list

    def deploy(self):
        #raw_bam = self.find_bam(self.rbamDir, "^pick.*bam$")
        duplex_bam = self.find_bam(self.cbamDir, "^pick.*uBAM$")
        duplex_bam.append(self.cbam)
        duplex_number = self.run_count(duplex_bam)
        #raw_reads_number = self.run_count(raw_bam)
        print(duplex_number)
        full_raw_reads = self.run_count([self.rbam])[0]
        full_raw_reads = int(full_raw_reads)
        raw_reads_number = [full_raw_reads]
        for i in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]:
            raw_reads_number.append(int(full_raw_reads * i))
        raw_reads_number.sort()
        print(raw_reads_number)
        #self.plot(raw_reads_number, duplex_number)
        lowess_value = lowess(raw_reads_number, duplex_number)
        print(lowess_value)
        duplex_1 = lowess_value.T[0][-1]
        duplex_2 = lowess_value.T[0][4]     #duplex 0.5
        PCTratio = round(duplex_2 / duplex_1, 3)
        #lowess拟合后的曲线
        self.plot(lowess_value.T[1], lowess_value.T[0], duplex_1, duplex_2)
        return PCTratio

class PCT_status:
    def __init__(self, ):
        pass

    def judge_status(self, PCTratio):
        status = ""
        if PCTratio >= 0.9:
            status = "Pass"
        elif 0.7 <= PCTratio < 0.9:
            status = "Warnning"
        else:
            status = "Fail"
        return status

    def write(self, PCTratio, status, outfile):
        out = open(outfile, "w")
        out.write("PCTratio\tStatus\n")
        out.write(str(PCTratio) + "\t" + status+"\n")
        out.close()

    def deploy(self, PCTratio, outfile):
        status = self.judge_status(PCTratio)
        self.write(PCTratio, status, outfile)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="UMI qc")
    parser.add_argument('--cbam', help="duplex consensus bam file", required=True)
    parser.add_argument('--rbam', help="bam for call duplex consensus reads", required=True)
    parser.add_argument('--rbamDir', help="rbam directory", required=True)
    parser.add_argument('--cbamDir', help="cbam directory", required=True)
    parser.add_argument('--outfile', help="output file", required=True)
    parser.add_argument('--png', help="output png file")
    parser.add_argument('--samtools', help="samtools", default="/mnt/share01/tools/bin/samtools")
    parser.add_argument('--cpu', help="cpu number", type=int, default=8)
    argv = parser.parse_args()
    #计算抽取bam的reads数和生成duplex的reads数
    t = lowess_estimate(argv.rbamDir, argv.cbamDir, argv.rbam, argv.cbam, argv.samtools, argv.cpu, argv.png)
    PCTratio = t.deploy()
    #计算ratio，并判断质控是否通过
    t = PCT_status()
    t.deploy(PCTratio, argv.outfile)
