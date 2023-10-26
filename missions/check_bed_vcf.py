#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
室间质评：bed与vcf匹配，如果在vcf中，修改bed并加一列check，否则加nocheck
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-12-10 10:14:57"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import glob
import gzip
import argparse
import pandas as pd
import numpy as np


def trans_float(x):
    if x != '-':
        x = format(float(x.strip('%')) / 100, '.2f')
    return x


def trans_type(x):
    """
    Variant_Type_cat列内容：
    missense,C,A
    移码突变,A,-
    """
    vartype, ref, alt = x.split(',')
    if vartype == '移码突变':
        if ref == '-':
            x = 'insertion'
        elif alt == '-':
            x = 'deletion'
        else:
            x = 'complex'
    else:
        x = vartype
    return x


def trans_check(x):
    if x != 'no_check':
        x = 'check'
    return x


def get_bed_library_df(file_path):
    df = pd.read_table(file_path, dtype={'Start_Position': 'str'})
    df['Position'] = df['Chromosome'].str.cat(df['Start_Position'], sep='_')
    df['End_Position'] = df['Start_Position']
    df['VAF'] = df['VAF'].apply(trans_float)
    df['Variant_Type'] = df['Variant_Type'].str.replace('同义突变', 'silent', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('错义突变', 'missense', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('非编码区突变', 'non-coding', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('非移码突变', 'complex', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('移码缺失', 'deletion', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('剪切位点突变', 'splice', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('截断突变', 'nonsense', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('延长突变', 'stoploss', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('剪切位点附近的突变', 'non-coding', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('起始密码子突变', 'nonsense', regex=False)
    df['Variant_Type'] = df['Variant_Type'].str.replace('终止密码子突变', 'stoploss', regex=False)
    df['Variant_Type_cat'] = df['Variant_Type'].str.cat([df['Reference_Allele'], df['Alternative_Allele']], sep=',')
    df['Variant_Type_cat'] = df['Variant_Type_cat'].apply(trans_type)
    df['Variant_Type'] = df['Variant_Type_cat']
    return df


def get_vcf_library_df(file_path):
    # tmp = f'{os.path.dirname(file_path)}/get_positive_mutation.tmp.txt'
    tmp = 'get_positive_mutation.tmp.txt'
    tmp_w = open(tmp, 'w')
    if 'gz' in file_path:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp_w.write(line)
    else:
        with open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp_w.write(line)
    tmp_w.close()
    df = pd.read_table(tmp, dtype={'POS': 'str'})
    os.remove(tmp)

    df['Position'] = df['#CHROM'].str.cat(df['POS'], sep='_')
    df = df[['Position', '#CHROM']]
    return df


def get_args():
    parser = argparse.ArgumentParser(description='bed与vcf匹配，如果在vcf中，修改bed并加一列check，否则加nocheck')
    parser.add_argument('--bed', help='bed文件路径', required=True)
    parser.add_argument('--vcf', help='vcf文件路径', required=True)
    parser.add_argument('--out', help='输出文件前缀', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    pd.set_option('max_columns', 15)
    argv = get_args()
    bed_df = get_bed_library_df(argv.bed)
    vcf_df = get_vcf_library_df(argv.vcf)
    merge_df = bed_df.merge(vcf_df, how='left', on='Position')
    merge_df['#CHROM'] = merge_df['#CHROM'].fillna('no_check')
    merge_df['#CHROM'] = merge_df['#CHROM'].apply(trans_check)
    # merge_df['Chromosome'] = merge_df['Chromosome'].str.strip('chr')
    merge_df = merge_df.rename(columns={'Chromosome': '##Chromosome'})
    merge_df = merge_df.drop(columns=['Position', 'Variant_Type_cat'])
    check_df = merge_df[merge_df['#CHROM'] == 'check']
    nocheck_df = merge_df[merge_df['#CHROM'] == 'no_check']
    merge_df = merge_df.drop(columns=['#CHROM'])
    check_df = check_df.drop(columns=['#CHROM'])
    nocheck_df = nocheck_df.drop(columns=['#CHROM'])
    merge_df.to_csv(f'{argv.out}.bed', sep='\t', index=False)
    check_df.to_csv(f'{argv.out}.check.xls', sep='\t', index=False)
    nocheck_df.to_csv(f'{argv.out}.no_check.xls', sep='\t', index=False)

