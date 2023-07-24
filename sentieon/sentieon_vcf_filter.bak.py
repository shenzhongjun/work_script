#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sentieon *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断
输入：*.filtered.tmp.vcf.gz路径、*.filtered.vcf.gz路径、vcf类型（双样本、单样本、白名单）
输出：vcf.gz格式过滤结果
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-4-1 13:48:11"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import re
import pandas as pd
import argparse
import gzip


def get_args():
    parser = argparse.ArgumentParser(description='Sentieon *.filtered.tmp.vcf.gz结果依据filter列标签进行可靠性判断。')
    parser.add_argument('--input_vcf', '-i', help='vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf.gz输出文件路径', required=True)
    parser.add_argument('--vcf_type ', '-t', help='vcf类型，配对样本、单样本或白名单等')
    return parser.parse_args()


def filter_vcf(x):
    pas = 'PASS'
    nopass = 'NOPASS'
    uncert = 'UNCERTAIN'

    qual = int(x['qual'])
    reads = int(x['tumor_id'].split(':')[1].split(',')[1])
    vaf = float(x['tumor_id'].split(':')[2])
    normal_vaf = float(x['normal_id'].split(':')[2])
    low_vcf_ratio = (vaf > 0 and normal_vaf > 0 and float(vaf/normal_vaf) < 3)
    noisy = int(re.search("ECNT=(.*?);", x['info']).group(1))
    if 'SOR=' in x['info']:
        strand = float(re.search("SOR=(.*?);", x['info']).group(1))
    else:
        strand = 0
    db = x['info'].startswith('DB')
    if len(x['tumor_id'].split(':')) >= 14:
        fusite = '|' in x['tumor_id'].split(':')[13]
    else:
        fusite = False

    # 过滤规则：
    def filter_strand_bias():
        if db and qual >= 100 and reads >= 8 and vaf >= 0.01 and noisy <= 5:
            return nopass  # uncert
        elif qual >= 500 and reads >= 8 and vaf >= 0.01 and noisy <= 5:
            return nopass  # uncert
        else:
            return nopass

    def filter_noisy_region():
        if not fusite and qual >= 200 and vaf >= 0.01 and reads >= 8:
            return nopass  # uncert
        else:
            return nopass

    def filter_short_tandem_repeat():
        if qual >= 500 and reads >= 8 and vaf >= 0.01 and noisy <= 5 and normal_vaf <= 0.001:
            return pas
        elif qual >= 100 and reads >= 8 and vaf >= 0.01 and noisy <= 10 and normal_vaf <= 0.001:
            return nopass  # uncert
        else:
            return nopass

    def filter_orientation_bias():
        if db and qual >= 10 and reads >= 8 and vaf >= 0.01 and noisy <= 5:
            return pas
        elif qual >= 100 and reads >= 8 and vaf >= 0.01 and noisy <= 5:
            return pas
        elif db and qual >= 1 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        elif qual >= 10 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        else:
            return nopass

    def filter_read_pos_bias():
        if db and qual >= 10 and reads >= 8 and vaf >= 0.01 and noisy <= 5 and not fusite:
            return pas
        elif qual >= 100 and reads >= 8 and vaf >= 0.01 and noisy <= 5 and not fusite:
            return pas
        elif db and qual >= 1 and reads >= 8 and vaf >= 0.005 and noisy <= 10 and not fusite:
            return nopass  # uncert
        elif qual >= 10 and reads >= 8 and vaf >= 0.005 and noisy <= 10 and not fusite:
            return nopass  # uncert
        else:
            return nopass

    def filter_low_qual():
        if qual >= 5 and reads >= 8 and vaf >= 0.01 and noisy <= 5 and strand <= 1:
            return pas
        elif qual >= 1 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        else:
            return nopass

    def filter_triallelic_site():
        if db and qual >= 10 and reads >= 8 and vaf >= 0.005 and noisy <= 5:
            return pas
        elif qual >= 100 and reads >= 8 and vaf >= 0.007 and noisy <= 5:
            return pas
        elif db and qual >= 1 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        elif qual >= 10 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        else:
            return nopass

    def filter_insignificant():     # 如果normal_vcf不为0，加入频率比
        if qual >= 10 and reads >= 8 and vaf >= 0.005 and noisy <= 5:
            return pas
        elif qual >= 1 and reads >= 8 and vaf >= 0.005 and noisy <= 10:
            return nopass  # uncert
        else:
            return nopass

    # 过滤标签按可靠性风险分为三类：风险tag（更容易被判断为假阳）、一般tag（有可能是假阳）和标志tag（不影响变异可靠性）
    filter_tags = set(x['filter'].split(';'))
    risk_tags = {'strand_bias', 'noisy_region', 'short_tandem_repeat'}
    normal_tags = {'orientation_bias', 'read_pos_bias', 'low_qual', 'triallelic_site', 'insignificant'}
    sign_tags = {'germline_risk', 'alt_allele_in_normal', 'low_t_alt_frac'}
    risk_in_filter_tags = risk_tags.intersection(filter_tags)
    normal_in_filter_tags = normal_tags.intersection(filter_tags)
    sign_in_filter_tags = sign_tags.intersection(filter_tags)

    # 如果是PASS则跳过
    if pas in filter_tags:
        pass
    # 如果不是PASS，首先判断是否融合断点
    elif fusite:
        filter_tags.add('fusion_clip')
        filter_tags.add(nopass)
    else:
        # 依据风险tag和一般tag进行过滤
        if len(risk_in_filter_tags) == 1:
            filter_tags.add(eval(f'filter_{list(risk_in_filter_tags)[0]}')())
        elif len(risk_in_filter_tags) >= 2:
            filter_tags.add(nopass)
        elif len(normal_in_filter_tags) >= 1:
            filter_result = []
            for tag in normal_in_filter_tags:
                filter_result.append(eval(f'filter_{tag}')())
            if nopass in filter_result:
                filter_tags.add(nopass)
            elif uncert in filter_result:
                filter_tags.add(uncert)
            else:
                filter_tags.add(pas)
        elif len(sign_in_filter_tags) >= 1:
            filter_tags.add(pas)
        else:
            pass

    # 处理标志tag，优化germline_risk判断
    if 't_lod_fstar' in filter_tags:
        filter_tags.remove('t_lod_fstar')
    if (('alt_allele_in_normal' in filter_tags) or ('insignificant' in filter_tags)) and low_vcf_ratio:
        filter_tags.add('germline_risk')
    x['filter'] = ';'.join(sorted(list(filter_tags)))
    return x


if __name__ == "__main__":
    args = get_args()
    input_vcf = args.input_vcf
    output_vcf = args.output_vcf
    with gzip.open(input_vcf, 'rt') as f, gzip.open(output_vcf, 'wt') as w:
        for line in f:
            if line.startswith('#'):
                w.write(line)
    df = pd.read_table(
        input_vcf,
        compression='gzip',
        comment='#',
        names=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'tumor_id', 'normal_id']
    )
    df = df.apply(filter_vcf, axis=1)
    df.to_csv(output_vcf, mode='a', header=False, sep='\t', index=False, compression='gzip')

