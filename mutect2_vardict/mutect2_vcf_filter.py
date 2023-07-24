#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mutect2 *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断
输入：*.filtered.tmp.vcf.gz路径、*.filtered.vcf.gz路径、vcf类型（双样本、单样本、白名单）
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
    parser = argparse.ArgumentParser(description='Mutect2 *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断。')
    parser.add_argument('--input_vcf', '-i', help='vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf.gz输出文件路径', required=True)
    # parser.add_argument('--vcf_type ', '-t', help='vcf类型，配对样本、单样本或白名单等')
    return parser.parse_args()


def tag_vcf(x, t, n):
    """
    Mutect2 vcf FILTER INFO FORMAT TUMOR列示例：
    # normal_artifact;slippage;weak_evidence
    # AS_FilterStatus=weak_evidence;AS_SB_TABLE=252,271|10,13;DP=580;ECNT=1;GERMQ=93;MBQ=20,20;MFRL=172,159;MMQ=60,60;MPOS=33;NALOD=-0.4575;NLOD=19.13;POPAF=6;ROQ=93;RPA=8,7;RU=C;STR;STRQ=1;TLOD=3.87
    # GT:AD:AF:DP:F1R2:F2R1:FAD:SB
    # 0/1:411,17:0.026:428:140,6:139,5:310,12:199,212,8,9	0/0:112,6:0.041:118:34,1:45,3:86,4:53,59,2,4
    """
    filter_tags = set(x['filter'].split(';'))

    if 'germline' in filter_tags:
        filter_tags.remove('germline')
        filter_tags.add('germline_risk')

    ref_reads, reads = [int(a) for a in x[t].split(':')[1].split(',')]
    total_reads = int(x[t].split(':')[3])
    vaf = reads / total_reads
    mbq = float(re.search("MBQ=(.*?);", x['info']).group(1).split(',')[0]) if 'MBQ=' in x['info'] else 30
    mpos = int(re.search("MPOS=(.*?);", x['info']).group(1).split(',')[0]) if 'MPOS=' in x['info'] else 30
    msirep = int(re.search("RPA=(.*?);", x['info']).group(1).split(',')[0]) if 'RPA=' in x['info'] else 0

    noraml_ref_reads, normal_reads = [int(a) for a in x[n].split(':')[1].split(',')]
    normal_total_reads = int(x[n].split(':')[3])
    normal_total_reads = 1 if normal_total_reads == 0 else normal_total_reads
    normal_vaf = normal_reads / normal_total_reads
    high_normal_vaf = (normal_vaf >= 0.02)
    low_vaf_ratio = (vaf > 0 and normal_vaf > 0 and float(vaf / normal_vaf) < 4)

    # 严格以肿瘤对照比来判定突变来源不确定
    if low_vaf_ratio:       # (('normal_artifact' in filter_tags) or high_normal_vaf or filter_tags == {'PASS'}) and :
        filter_tags.add('germline_risk')

    # #### 过滤 ####
    # 过滤低于检出限的突变
    if vaf < 0.006 or reads < 4:        # Mutect2实际上无力检出0.5%频率的突变，一般检出限在1%，甚至个别2%的也会漏检
        filter_tags.add('below_lod')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # 过滤位置偏好性突变
    if mpos <= 12 and (vaf < 0.01 or normal_reads > 1):
        filter_tags.add('position_bias')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # 过滤同一单体型上的突变
    if 'PGT' in x['format'] and (vaf < 0.01 or normal_reads > 1):
        filter_tags.add('haplotype')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # 严格过滤STR
    if (msirep >= 3 and (vaf <= 0.05 or normal_reads >= 1)) or (msirep >= 5 and (vaf <= 0.1 or normal_reads >= 1 or normal_total_reads <= 100)):
        filter_tags.add('short_tandem_repeat')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # 过滤低频且质量不高的indel
    if (len(x['ref']) >= 2 or len(x['alt']) >= 2) and vaf <= 0.01 and mbq < 30:
        filter_tags.add('low_mbq')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # 大于阈值50bp的indel不报出
    if len(x['ref']) > 50 or len(x['alt']) > 50:
        filter_tags.add('long_indel')
        if 'PASS' in filter_tags:
            filter_tags.remove('PASS')

    # #### 捞回 ####
    can_pass_tags = {'multiallelic', 'slippage', 'normal_artifact', 'germline_risk'}
    if filter_tags.issubset(can_pass_tags):
        # Mutect2检出的多等位基因突变多位于STR区域或旁边，indel不可信，与normal_artifact同时出现也不可信
        if 'multiallelic' in filter_tags:
            if 'normal_artifact' not in filter_tags and len(x['ref']) == 1 and len(x['alt']) == 1:
                filter_tags.add('PASS')
        else:
            filter_tags.add('PASS')

    x['filter'] = ';'.join(sorted(list(filter_tags)))
    return x


def check_header(header, t, n):
    cols = header.split('\t')[:9]
    sample1, sample2 = header.split('\t')[9:]

    if sample1 == t and sample2 == n:
        return [header, False]
    elif sample1 == n and sample2 == t:
        print(f'======== VCF中肿瘤和对照样本位置与预期顺序不符，已修改为肿瘤在前对照在后 ======== ')
        return ['\t'.join(cols + [t, n]), True]
    else:
        print(f'tumor="{t}", normal={n}. sample1 in vcf={sample1}，sample2 in vcf={sample2}')
        raise BaseException(f'======== VCF的样本与给出的样本名不符，请核对！ ======== ')


if __name__ == "__main__":
    args = get_args()
    input_vcf = args.input_vcf
    output_gzvcf = args.output_vcf
    output_vcf = f'{output_gzvcf.split(".gz")[0]}'
    change_headerline = False
    tumor = ''
    normal = ''

    with gzip.open(input_vcf, 'rt') as f, open(output_vcf, 'wt') as w:
        for line in f:
            if line.startswith('##FILTER=<ID=germline'):
                pass
            elif line.startswith('##tumor_sample'):
                tumor = line.strip().split('=')[1]
                w.write(line)
            elif line.startswith('##normal_sample'):
                normal = line.strip().split('=')[1]
                w.write(line)
            elif line.startswith('##'):
                w.write(line)
            elif line.startswith('#CHROM'):
                header = line.strip()

        w.write('##FILTER=<ID=germline_risk,Description="Evidence indicates this site is germline, not somatic">\n')
        w.write('##FILTER=<ID=below_lod,Description="Far below limit of detection">\n')
        # w.write('##FILTER=<ID=indel_site,Description="Variation is next to an indel">\n')
        w.write('##FILTER=<ID=position_bias,Description="Variation position is close to reads end">\n')
        w.write('##FILTER=<ID=low_mbq,Description="Low mbq indel with low af">\n')
        w.write('##FILTER=<ID=short_tandem_repeat,Description="Short tandem repeat, may be caused by sequencing error">\n')
        w.write('##FILTER=<ID=long_indel,Description="Long Indel more than 50bp, which difficult to verify">\n')
        w.write(f'##Mutect2_vcf_filterCommand="python {os.path.abspath(__file__)} -i {input_vcf} -o {output_gzvcf}"\n')
        header, change_headerline = check_header(header, tumor, normal)
        w.write(f'{header}\n')

    if change_headerline:
        names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', normal, tumor]
    else:
        names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', tumor, normal]
    df = pd.read_table(
        input_vcf,
        compression='gzip',
        comment='#',
        names=names
    )
    if change_headerline:
        df = df.reindex(
            ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', tumor, normal], axis=1
        )
    df = df.apply(tag_vcf, axis=1, args=(tumor, normal))
    df.to_csv(output_vcf, mode='a', header=False, sep='\t', index=False)
    pysam.tabix_compress(output_vcf, output_gzvcf, force=True)
    pysam.tabix_index(output_gzvcf, force=True, preset='vcf')
    os.remove(output_vcf)
