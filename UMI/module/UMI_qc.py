# /usr/bin/env python
# conding:utf-8
import os
import subprocess
import argparse
import pysam
from multiprocessing import Pool

def read_bed(bed):
    beddic = []
    for line in open(bed, "r").readlines():
        beddic.append(line.strip().split("\t"))
    return beddic

def QC(bam, chrom, start, end):
    bamFile = pysam.AlignmentFile(bam, "rb")
    UMI_index = {}
    for read in bamFile.fetch(chrom, int(start), int(end)):
        tag = read.get_tag("RX")
        if "N" in tag:
            continue
        left, right = tag.split("-")
        if read.is_read1:
            UMI_index.setdefault(left, 0)
            UMI_index[left] += 1
        else:
            UMI_index.setdefault(right, 0)
            UMI_index[right] += 1
    bamFile.close()
    #print(UMI_index)
    return UMI_index

def write_output(dict, outfile):
    out = open(outfile, "w")
    """
    #可以根据需要修改输出的文件格式
    out.write("UMItype\tUMN_Number\n")
    for key in sorted(dict.keys()):
        out.write(key + "\t" + str(dict[key]) + "\n")
    """
    out.write("分子标签统计\t统计结果\n")
    out.write("Total_Type\t" + str(len(dict.keys())) + "\n")
    max_value = max(dict.values())
    min_value = min(dict.values())
    Fold = round(max_value / min_value, 2)
    QC_status = ""
    if Fold < 8:
        QC_status = "Pass"
    elif 8 <= Fold < 20:
        QC_status = "Warnning"
    elif Fold >= 20:
        QC_status = "Fail"
    out.write("max_count\t" + str(max_value) + "\n")
    out.write("min_count\t" + str(min_value) + "\n")
    out.write("Fold\t" + str(Fold) + "\n" )
    out.write("QC_result\t" + QC_status + "\n")
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="UMI qc")
    parser.add_argument('--input', help="input bam file", required=True)
    parser.add_argument('--output', help="output file", required=True)
    parser.add_argument('--cpu', help="cpu number", type=int, default=8)
    parser.add_argument('--bed', help="bed")
    argv = parser.parse_args()

    pool = Pool(argv.cpu)
    result = []
    beddic = read_bed(argv.bed)
    for info in beddic:
        chrom, start, end = info
        result.append(pool.apply_async(QC, (argv.input, chrom, start, end, )))
    pool.close()
    pool.join()

    UMI_index = {}
    for res in result:
        if not UMI_index:
            UMI_index = res.get()
        else:
            for key, value in res.get().items():
                UMI_index.setdefault(key, 0)
                UMI_index[key] += value
    #print(UMI_index)
    write_output(UMI_index, argv.output)



