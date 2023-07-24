#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate TMB value by varscan + mutect2 annotation files.
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-9-8 14:02:20"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os

import numpy as np
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='Calculate TMB value by varscan + mutect2 annotation files.')
    parser.add_argument('--comm', help='the varscan & mutect2 comm vcf file.', required=True)
    parser.add_argument('--uniq', help='the varscan uniq vcf file.', required=True)
    return parser.parse_args()


def filter_by_keyword(obj, column, word):
    return obj[column].str.contains(word) is True


def filter_by_spv(obj, column, spv):
    return obj[column].str.rsplit(';', expand=True, n=1)[1].str.split('=', expand=True)[1].astype(np.float64) <= spv


def filter_by_depth(obj, sample, depth):
    return obj[sample].str.split(':', expand=True)[2].astype(np.int64) > depth


def filter_by_vardepth(obj, sample, depth):
    if sample == 'TUMOR':
        return obj[sample].str.split(':', expand=True)[4].astype(np.int64) >= depth
    elif sample == 'NORMAL':
        return obj[sample].str.split(':', expand=True)[4].astype(np.int64) < depth


def filter_by_freq(obj, sample, freq):
    if sample == 'TUMOR':
        return obj[sample].str.split(':', expand=True)[5].str.split('%', expand=True)[0].astype(np.float64) >= freq
    elif sample == 'NORMAL':
        return obj[sample].str.split(':', expand=True)[5].str.split('%', expand=True)[0].astype(np.float64) < freq


def filter_var_mu2_comm(obj):
    """
    related head lines in vcf file:
    ##INFO=<ID=SS,Number=1,Type=String,Description="Somatic status of variant (0=Reference,1=Germline,2=Somatic,3=LOH, or 5=Unknow)">
    ##INFO=<ID=SPV,Number=1,Type=Float,Description="Fisher's Exact Test P-value of tumor versus normal for Somatic/LOH calls">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
    ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Depth of variant-supporting bases (reads2)">
    ##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">

    the 'INFO' is like:
    DP=1917;SOMATIC;SS=2;SSC=15;GPV=1E0;SPV=2.5275E-2

    the 'FORMAT' is like:
    GT:GQ:DP:RD:AD:FREQ:DP4

    the 'TUMOR' is like:
    0/0:.:971:962:9:0.93%:627,335,7,2
    """
    df = obj[(filter_by_keyword(obj, 'INFO', 'SS=2')) & (filter_by_spv(obj, 'INFO', 0.01))]
    df = df[(filter_by_depth(df, 'TUMOR', 30)) & filter_by_depth(df, 'NORMAL', 10)]
    df = df[(filter_by_vardepth(df, 'TUMOR', 5)) & (filter_by_vardepth(df, 'NORMAL', 5))]
    df = df[(filter_by_freq(df, 'TUMOR', 1)) & (filter_by_freq(df, 'NORMAL', 1))]
    return df


def filter_varscan_uniq(obj):
    df = obj[filter_by_spv(obj, 'INFO', 0.01)]
    df = df[(filter_by_freq(df, 'TUMOR', 5)) & (-filter_by_freq(df, 'NORMAL', 5))]
    return df


if __name__ == "__main__":
    pd.set_option('max_rows', 5)
    pd.set_option('max_columns', 15)
    args = get_args()
    comm = pd.read_table(args.comm, header=17)
    uniq = pd.read_table(args.uniq, header=17)
    tmb_num = len(comm) + len(uniq)
    tmb_value = f'{tmb_num / 39:.2f}'
    path = os.path.split(os.path.realpath(args.anno))[0]
    stat = pd.DataFrame({'TMB_somatic_count': [tmb_num], 'TMB': [tmb_value]})
    stat_name = f'{os.path.splitext(os.path.basename(args.anno))[0]}.TMB_stat.xls'
    stat.to_csv(os.path.join(path, stat_name), sep='\t', index=False)
    print(stat)
    print(filter_var_mu2_comm(comm).shape)
    print(filter_varscan_uniq(uniq).shape)
