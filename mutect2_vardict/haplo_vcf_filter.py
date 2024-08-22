#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
GATK Haplotype *.haplotype.raw.vcf.gz进行基础过滤
输入：*.haplotype.raw.vcf.gz
输出：vcf.gz格式过滤结果
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-11-7 8:48:11"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import re
import pandas as pd
import argparse
import gzip
import pysam


def get_args():
    parser = argparse.ArgumentParser(description='GATK Haplotype *.haplotype.raw.vcf.gz进行基础过滤。')
    parser.add_argument('--input_vcf', '-i', help='vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf.gz输出文件路径', required=True)
    parser.add_argument('--tumor', '-t', help='肿瘤样本名', required=True)
    parser.add_argument('--normal', '-n', help='对照样本名', required=True)
    parser.add_argument('--tumor2', '-t2', help='WES+Panel检测中的Panel肿瘤样本名')
    return parser.parse_args()


def tag_vcf(x, t, n):
    """
    Haplotype vcf FILTER INFO FORMAT AJL105 AJL105DZ AJL106(WES_T,WES_N,PANEL_T)列示例：
    # .
    # AC=6;AF=1.00;AN=6;DP=276;ExcessHet=0.0000;FS=0.000;MLEAC=5;MLEAF=0.833;MQ=60.00;QD=25.36;SOR=0.959
    # GT:AD:DP:GQ:PL
    # 1/1:0,233:233:99:10332,701,0  1/1:0,34:34:99:1530,102,0   1/1:0,1:1:3:45,3,0
    # 原bcftools过滤条件'REF="*" || ALT="*" || SUM(ALT)>=30 || SUM(REF)>=30 || QD<1 || MQ<20 || FS>100 || SOR>20 ||
    ReadPosRankSum<-20 || FMT/DP<20 || FMT/AD[0:1]<5 || (FMT/AD[0:1])/(FMT/DP)<0.25'
    """
    x['filter'] = 'PASS'
    if './.' in x[t] or './.' in x[n]:
        x['filter'] = '.'
    else:
        ref_reads, reads = [int(a) for a in x[t].split(':')[1].split(',')]
        total_reads = int(x[t].split(':')[2])
        total_reads = 1 if total_reads == 0 else total_reads
        vaf = reads / total_reads

        noraml_ref_reads, normal_reads = [int(a) for a in x[n].split(':')[1].split(',')]
        normal_total_reads = int(x[n].split(':')[2])
        normal_total_reads = 1 if normal_total_reads == 0 else normal_total_reads
        normal_vaf = normal_reads / normal_total_reads

        fs = float(re.search("FS=(.*?);", x['info']).group(1)) if 'FS=' in x['info'] else 0
        mq = float(re.search("MQ=(.*?);", x['info']).group(1)) if 'MQ=' in x['info'] else 60
        qd = float(re.search("QD=(.*?);", x['info']).group(1)) if 'QD=' in x['info'] else 0
        rpos = float(re.search("ReadPosRankSum=(.*?);", x['info']).group(1)) if 'ReadPosRankSum=' in x['info'] else 0
        sor = float(re.search("SOR=(.*)", x['info']).group(1)) if 'SOR=' in x['info'] else 0

        if (x['ref'] == '*' or x['alt'] == '*' or len(x['ref']) >= 30 or len(x['alt']) >= 30 or fs > 80 or mq < 40 or
                qd < 1 or rpos < -20 or sor > 20 or normal_vaf < 0.25 or normal_total_reads < 15 or normal_reads < 5 or
                total_reads < 20 or reads < 5):
            x['filter'] = '.'
    return x


if __name__ == "__main__":
    args = get_args()
    input_vcf = args.input_vcf
    output_gzvcf = args.output_vcf
    output_vcf = f'{output_gzvcf.split(".gz")[0]}'
    tumor = args.tumor
    normal = args.normal
    tumor2 = args.tumor2

    with gzip.open(input_vcf, 'rt') as f, open(output_vcf, 'wt') as w:
        for line in f:
            if line.startswith('##'):
                w.write(line)
            elif line.startswith('#CHROM'):
                header = line
                vcf_samples = line.strip().split('\t')[9:]
        w.write(f'##tumor_sample={tumor}\n')
        if tumor2:
            w.write(f'##tumor_sample2={tumor2}\n')
        w.write(f'##normal_sample={normal}\n')
        w.write(f'##Haplotype_vcf_filterCommand="python {os.path.abspath(__file__)} -i {input_vcf} -o {output_gzvcf} -t {tumor} -n {normal}"\n')
        w.write(header)

    names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + vcf_samples
    df = pd.read_table(
        input_vcf,
        compression='gzip',
        comment='#',
        names=names
    )

    df = df.apply(tag_vcf, axis=1, args=(tumor, normal))
    df = df.loc[df['filter'] == 'PASS']
    df.to_csv(output_vcf, mode='a', header=False, sep='\t', index=False)
    pysam.tabix_compress(output_vcf, output_gzvcf, force=True)
    pysam.tabix_index(output_gzvcf, force=True, preset='vcf')
    os.remove(output_vcf)

