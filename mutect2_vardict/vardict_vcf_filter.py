#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Vardict *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断
输入：*.filtered.tmp.vcf.gz路径
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
    parser = argparse.ArgumentParser(description='Vardict *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断。')
    parser.add_argument('--input_vcf', '-i', help='vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf.gz输出文件路径', required=True)
    return parser.parse_args()


def tag_vcf(x):
    """
    Vardict vcf FILTER INFO FORMAT TUMOR列示例：
    # P0.05
    # STATUS=StrongSomatic;SAMPLE=AHT653;TYPE=SNV;DP=932;VD=4;AF=0.0043;SHIFT3=0;MSI=1;MSILEN=1;SSF=0.34231;SOR=0;LSEQ=CTACAATTATGTGTGGCAAC;RSEQ=CATGTGTCGGCATTATGGCT
    # GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM
    # 0/1:932:4:2,2:491,437:928,4:0.0043:2,2:38.2:1:37.5:1:1:1.12343:58:8:0.0043:0.0011:2
    """
    filter_tags = set(x['filter'].split(';'))
    infos = x['info'].split(';')
    qual = int(x['qual'])
    status = infos[0].split('=')[1]
    depth = int(infos[3].split('=')[1])
    reads = int(infos[4].split('=')[1])
    vaf = float(infos[5].split('=')[1])
    msirep = float(infos[7].split('=')[1])
    msilen = int(infos[8].split('=')[1])
    msi = (msirep >= 3)
    indellikely = 'InDelLikely' in x['info']
    p001likely = 'P0.01Likely' in x['info']
    mut_type = infos[2].split('=')[1]
    normal_reads = int(x['normal_id'].split(':')[2])
    normal_vaf = float(x['normal_id'].split(':')[6])
    high_normal_vaf = (normal_vaf >= 0.02)
    low_vcf_ratio = (vaf > 0 and normal_vaf > 0 and float(vaf / normal_vaf) < 4)

    # 对P0.05的SNV变异进行捞回
    if mut_type == 'SNV' and vaf >= 0.004 and reads >= 4 and normal_reads == 0 and not msi and filter_tags == {'P0.05'}:
        filter_tags.add('PASS')

    # 过滤低于检出限的突变
    if vaf < 0.004 or reads < 4:
        filter_tags.add('below_lod')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    if reads > depth and mut_type != 'SNV':     # 处理突变深度大于总深度的情况
        filter_tags.add('complex_indel')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    if len(x['ref']) > 30 or len(x['alt']) > 30:        # Vardict大于30bp的indel不报出
        filter_tags.add('long_indel')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    if status == 'StrongSomatic':
        if mut_type != 'SNV' and vaf >= 0.9:        # 处理超高频体细胞假阳性indel
            filter_tags.add('germline_risk')
            if 'PASS' in filter_tags:
                filter_tags.remove('PASS')
        elif reads <= 4 and x['format'] == 'GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM' \
                and float(x['tumor_id'].split(':')[-9]) <= 30:         # 处理支持reads少且质量不高的变异
            if 'PASS' in filter_tags:
                filter_tags.remove('PASS')
            filter_tags.add('low_qual')

    if status == 'Germline' or status == 'LikelyLOH':
        if normal_vaf >= 0.25:  # 软件判断为Germline或LikelyLOH基础上，对照突变频率25%以上判断为胚系变异，以下判断为突变来源不确定
            if filter_tags == {'P0.05'}:  # and not msi:
                filter_tags.add('PASS')
            filter_tags.add('germline')
        elif low_vcf_ratio:
            filter_tags.add('germline_risk')
    elif ((status == 'LikelySomatic' or high_normal_vaf or indellikely) and low_vcf_ratio) or (p001likely and (vaf <= 0.01 or normal_reads > 0)):
        filter_tags.add('germline_risk')

    if mut_type != 'SNV' and msi and 'germline' not in filter_tags and (vaf <= 0.05 or qual <= 100 or reads <= 10 or normal_reads > 0 or (msirep >= 5 and (vaf <= 0.1 or normal_reads > 0))):
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')
        filter_tags.add('short_tandem_repeat')

    x['filter'] = ';'.join(sorted(list(filter_tags)))
    return x


if __name__ == "__main__":
    args = get_args()
    input_vcf = args.input_vcf
    output_gzvcf = args.output_vcf
    output_vcf = f'{output_gzvcf.split(".gz")[0]}'

    with gzip.open(input_vcf, 'rt') as f, open(output_vcf, 'wt') as w:
        for line in f:
            if line.startswith('##'):
                w.write(line)
            if line.startswith('#CHROM'):
                header = line

        w.write('##FILTER=<ID=low_qual,Description="Low mean quality score in variant depth">\n')
        w.write('##FILTER=<ID=below_lod,Description="Far below limit of detection">\n')
        w.write('##FILTER=<ID=complex_indel,Description="Complex Indel with VD>DP, will change to equal manually">\n')
        w.write('##FILTER=<ID=long_indel,Description="Long Indel more than 50bp, which difficult to verify">\n')
        w.write('##FILTER=<ID=short_tandem_repeat,Description="Short tandem repeat, may be caused by sequencing error">\n')
        w.write('##FILTER=<ID=germline,Description="Evidence show this site is germline">\n')
        w.write('##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">\n')
        w.write(f'##Vardict_vcf_filterCommand="python {os.path.abspath(__file__)} -i {input_vcf} -o {output_gzvcf}"\n')
        w.write(header)

    df = pd.read_table(
        input_vcf,
        compression='gzip',
        comment='#',
        names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'tumor_id', 'normal_id']
    )
    df = df.apply(tag_vcf, axis=1)
    df.to_csv(output_vcf, mode='a', header=False, sep='\t', index=False)
    pysam.tabix_compress(output_vcf, output_gzvcf, force=True)
    pysam.tabix_index(output_gzvcf, force=True, preset='vcf')
    os.remove(output_vcf)
