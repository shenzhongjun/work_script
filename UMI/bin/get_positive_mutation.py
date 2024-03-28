#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
从UMI注释结果中提取阳性位点Call变异情况
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-10-29 10:16:38"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import glob
import re

import pandas as pd
import argparse
# pd.set_option('max_columns', 15)


def get_args():
    parser = argparse.ArgumentParser(description='从UMI注释结果中提取阳性位点Call变异情况。')
    parser.add_argument('--positive_list', help='阳性位点列表，要求包含"Position"列，内容为chr1_22222', required=True)
    parser.add_argument('--anno_path', help='UMI注释结果路径，支持"*"')
    parser.add_argument('--sen_vcf', help='Sentieon软件call变异结果路径，支持"*"')
    return parser.parse_args()


def get_file_paths(pathname):
    paths = glob.glob(pathname)
    p_dic = {}
    for full_path in paths:
        library_name = os.path.split(full_path)[1].split('.')[0]
        p_dic[library_name] = full_path
    return p_dic


def get_library_df(file_path, library_name):
    df = pd.read_table(file_path, dtype={'Start': 'str'})
    df['Position'] = df['Chr'].str.cat(df['Start'], sep='_')
    df = df[['Position', 'TotalDepth', 'RefDepth', 'AltDepth', 'MutFreq']]
    df = df.rename(columns={'TotalDepth': f'TotalDepth_{library_name}'})
    return df


def get_sen_paths(pathname):
    paths = glob.glob(pathname)
    p_dic = {}
    for full_path in paths:
        library_name = re.search(r'/DYDFZ-.*?/', full_path).group(0).split('/')[1]
        p_dic[library_name] = full_path
    return p_dic


def get_sen_df(file_path, library_name):
    df = pd.read_table(file_path, header=99, dtype={'POS': 'str'})
    df['Position'] = df['#CHROM'].str.cat(df['POS'], sep='_')
    df['RefDepth'] = df[library_name].str.split(':', expand=True)[1].str.split(',', expand=True)[0]
    df['AltDepth'] = df[library_name].str.split(':', expand=True)[1].str.split(',', expand=True)[1]
    df['MutFreq'] = df[library_name].str.split(':', expand=True)[2]
    df = df[['Position', 'RefDepth', 'AltDepth', 'MutFreq']]
    df = df.rename(columns={'RefDepth': f'RefDepth_{library_name}'})
    return df


if __name__ == "__main__":
    args = get_args()
    positive_list = args.positive_list
    mut_df = pd.read_table(os.path.join(positive_list))
    # anno_path example:
    # '/mnt/share05/clinical_project/projects/blood_tumor/other/zhufu_test/UMI/ctDNA*/DYDFZ*/step10.annotation/*hg19_multianno.reform.xls'
    if args.anno_path:
        path = args.anno_path
        mut_paths = get_file_paths(path)
        for library in mut_paths:
            library_df = get_library_df(mut_paths[library], library)
            mut_df = mut_df.merge(library_df, how='left', on='Position')
        mut_df = mut_df.drop(columns=['Position'])
        mut_df.to_csv(os.path.join(os.path.split(positive_list)[0], 'umi_positive_mutations.txt'), sep='\t', index=False)
    elif args.sen_vcf:
        path = args.sen_vcf
        mut_paths = get_sen_paths(path)
        for library in mut_paths:
            if os.path.exists(mut_paths[library]):
                library_df = get_sen_df(mut_paths[library], library)
                mut_df = mut_df.merge(library_df, how='left', on='Position')
        mut_df = mut_df.drop(columns=['Position'])
        mut_df.to_csv(os.path.join(os.path.split(positive_list)[0], 'sen_umi_positive_mutations.txt'), sep='\t', index=False)
