# /usr/bin/env python
# conding:utf-8
# 靶向捕获数据的质控
from argparse import ArgumentParser
import sys
import json
import pysam
import numpy as np
import subprocess
from multiprocessing import Pool

class target_QC:
    def __init__(self, cpu, bed, bam, fastp_json, samtools, outfile):
        self.cpu = cpu
        self.bed = bed
        self.bam = bam
        self.fastp_json = fastp_json
        self.samtools = samtools
        self.outfile = outfile

    def mapped_stat(self):
        """
        统计总的mapped reads和mapped bases
        """
        status, mapped_base = subprocess.getstatusoutput(fr"""{self.samtools} stats -s {self.bam}| grep 'bases mapped:' |cut -f 3 """)
        if status:
            print("统计mapping bases错误")
            sys.exit(1)
        status, dup_base = subprocess.getstatusoutput(fr"""{self.samtools} stats -s {self.bam} | grep 'bases duplicated:' |cut -f 3 """)
        if status:
            print("统计mapping bases错误")
            sys.exit(1)
        status, mapped_reads = subprocess.getstatusoutput(fr"""{self.samtools} flagstat {self.bam} | grep 'mapped' |grep "%"|cut -d " " -f 1 """)
        if status:
            print("统计mapping reads错误")
            sys.exit(1)
        status, total_reads = subprocess.getstatusoutput(fr"""{self.samtools} flagstat {self.bam} | grep 'total' |cut -d " " -f 1 """)
        if status:
            print("统计mapping reads错误")
            sys.exit(1)
        return mapped_base, mapped_reads, dup_base, total_reads

    def bam_stat(self, line):
        bamFile = pysam.AlignmentFile(self.bam, "rb")
        lines = line.split("\t")
        chrom, start, end = str(lines[0]), int(lines[1]), int(lines[2])
        bed_size = end - start
        ACGT = bamFile.count_coverage(chrom, start, end)
        target_bases = np.sum(ACGT)
        base_coverage_array = np.sum(np.array(ACGT), axis=0).tolist()
        ###获取PCR重复reads数
        dup_reads = 0
        target_reads = 0
        for alignment in bamFile.fetch(chrom, start, end, until_eof=True):
            if alignment.is_unmapped:
                continue
            target_reads += 1
            if alignment.is_duplicate:
                dup_reads += 1
        bamFile.close()
        return target_reads, target_bases, dup_reads, bed_size, base_coverage_array


    def read_json(self):
        infile = open(self.fastp_json, "r")
        json_info = json.load(infile)
        raw_bases = json_info['summary']['before_filtering']['total_bases']
        q20_rate = round(json_info['summary']['before_filtering']['q20_rate'] * 100, 2)
        q30_rate = round(json_info['summary']['before_filtering']['q30_rate'] * 100, 2)
        gc_content = round(json_info['summary']['before_filtering']['gc_content'] * 100, 2)
        dup_rate = round(json_info["duplication"]["rate"] * 100, 2)
        peak_insert_size = json_info["insert_size"]["peak"]
        adapter_bases = json_info["adapter_cutting"]["adapter_trimmed_bases"]
        adapter_rate = round(float(adapter_bases / raw_bases) * 100, 2)

        return [raw_bases, q20_rate, q30_rate, gc_content, dup_rate, peak_insert_size, adapter_rate]

    def deploy(self):
        raw_bases, q20_rate, q30_rate, gc_content, dup_rate, peak_insert_size, adapter_rate = self.read_json()
        mapped_base, mapped_reads, dup_base, raw_reads = self.mapped_stat()
        pool = Pool(self.cpu)
        result = []
        for line in open(self.bed, "r").readlines():
            line = line.strip()
            result.append(pool.apply_async(self.bam_stat, (line, )))
        pool.close()
        pool.join()
        target_reads, target_bases, dup_reads, bed_size = 0, 0, 0, 0
        base_coverage_array = []
        for res in result:
            t_read, t_base, dup_read, bed_s, base_coverage = res.get()
            target_reads += t_read
            target_bases += t_base
            dup_reads += dup_read
            bed_size += bed_s
            base_coverage_array.extend(base_coverage)
        #计算结果
        mapped_reads_ratio = round(float(int(mapped_reads) / int(raw_reads)) * 100, 2)
        mapped_data_ratio = round(float(int(mapped_base) / raw_bases) * 100, 2)
        PCR_dup_reads_ratio = round(float(dup_reads / int(raw_reads)) * 100, 2)
        target_data_ratio = round(float(target_bases / int(mapped_base)) * 100, 2)
        target_read_ratio = round(float(target_reads / int(mapped_reads)) * 100, 2)
        average_depth = round(target_bases / bed_size, 2)
        raw_depth = round(float(raw_bases / bed_size), 2)
        depth1 = average_depth * 0.2  # 0.2x平均深度
        depth2 = average_depth * 0.5  # 0.5x平均深度
        coverage_2 = np.round(float(np.sum([1 if i > depth1 else 0 for i in base_coverage_array]) / len(base_coverage_array)), 4)  # 0.2X平均深度的coverage
        coverage_5 = np.round(float(np.sum([1 if i > depth2 else 0 for i in base_coverage_array]) / len(base_coverage_array)), 4)  # 0.5X平均深度的coverage
        out = open(self.outfile, "w")
        out.write("Data_Size\tRaw_Reads\tRaw_Depth\tQ20(%)\tQ30(%)\tGC(%)\tDuplication(%)\tPeak_Insert_Size\tAdapter_Rate(%)\tMapped_Reads(%)\tMapped_Data(%)\tPCR_Duplicate_Reads(%)\tTarget_Data_in_Mapped_Data(%)\tTarget_Reads_in_Mapped_Reads(%)\tAverage_Depth(rmdup)\t0.2×_coverage\t0.5×_coverage\n")
        result = [str(int(raw_bases / 1000000))+"Mbp", raw_reads, raw_depth, q20_rate, q30_rate, gc_content, dup_rate, peak_insert_size, adapter_rate, mapped_reads_ratio, mapped_data_ratio, PCR_dup_reads_ratio, target_data_ratio, target_read_ratio, average_depth, coverage_2, coverage_5]
        result = [str(i) for i in result]
        out.write("\t".join(result)+"\n")

def main():
    parser = ArgumentParser()
    parser.add_argument("--sample", action="store", required=True, help="sample id")
    parser.add_argument("--bam", action="store", required=True, help="bam after markdup")
    parser.add_argument("--bed", action="store", required=True, help="regionfile or bed file")
    parser.add_argument("--json", action="store", required=True, help="json file of fastp")
    parser.add_argument("--qcFile", action="store", required=True, help="output qc file")
    parser.add_argument('--samtools', help="samtools", default="/mnt/share01/tools/bin/samtools")
    parser.add_argument("--cpu", action="store", default=8, type=int, help="number of cpus")
    argv = parser.parse_args()
    t = target_QC(argv.cpu, argv.bed, argv.bam, argv.json, argv.samtools, argv.qcFile)
    t.deploy()
if __name__ == '__main__':
    main()
