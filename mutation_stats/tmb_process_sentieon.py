#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量计算Sentieon vcf注释结果的TMB值，需要给出以SAMPLENAME为开头的文件路径。
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-9-7 10:24:41 "
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import glob
import argparse
import pandas as pd
import numpy as np
from io import StringIO


def get_file_paths(pathname):
    paths = glob.glob(pathname)
    p_dic = {}
    for full_path in paths:
        tumor = os.path.split(full_path)[1].split('.')[0]
        p_dic[tumor] = full_path
    return p_dic


def get_tumor_df(file_path):
    tmp2 = StringIO()
    with open(file_path, 'rt') as f:
        for line in f:
            if not line.startswith('##'):
                tmp2.write(line)
    df = pd.read_table(StringIO(tmp2.getvalue()), dtype={'POS': 'str'})
    return df


def filter_by_info(df):
    """根据vcf文件INFO列过滤"""
    # df = df[df['INFO'].str.contains('SOMATIC;SS=2')]    # filter SOMATIC
    filtered_df = df[df['INFO'].str.contains('nonsynonymous|stopgain|stoploss|nonframeshift')]
    # print(filtered_df)
    return filtered_df


def filter_by_format(df, tumor):
    """
    根据vcf文件FORMAT列过滤
    """
    # df_n = df['NORMAL'].str.split(':', expand=True)
    df_t = df[f'{tumor}'].str.split(':', expand=True).iloc[:, 1:3]
    # print(df_t)
    # df_n.columns = ['GT', 'GQ', 'DP', 'RD', 'AD', 'AF', 'DP4']
    df_t.columns = ['DP', 'AF']
    # df_n = df_n[(df_n['DP'].astype(np.int64) > 50) &
    #             (df_n['DP4'].str.split(',', expand=True)[2].astype(np.int64) +
    #              df_n['DP4'].str.split(',', expand=True)[3].astype(np.int64) >= 7)]
    # df_n = df_n[(df_n['DP'].astype(np.int64) > 50) &
    #             (df_n['AF'].str.split('%', expand=True)[0].astype(np.float64) <= 0.35)]
    # df_n = df_n[(df_n['DP'].astype(np.int64) > 50)]
    # df_t = df_t[(df_t['DP'].astype(np.int64) > 50) &
    #             (df_t['AF'].str.split('%', expand=True)[0].astype(np.float64) >= 5) &
    #             (df_t['DP4'].str.split(',', expand=True)[2].astype(np.int64) +
    #              df_t['DP4'].str.split(',', expand=True)[3].astype(np.int64) >= 8)]
    df_t = df_t[df_t['AF'].astype(np.float64) >= 0.05]
    # df_n_t = pd.concat([df_n, df_t], axis=1, join='inner')
    # df_final = df.loc[df_n_t.index.tolist()]
    df_final = df.loc[df_t.index.tolist()]
    return df_final


def calculate_tmb(glob_path, chip, out):
    mut_paths = get_file_paths(glob_path)
    result_dic = {}
    for tumor in mut_paths:
        df = get_tumor_df(mut_paths[tumor])
        filter_info = filter_by_info(df)
        filter_format = filter_by_format(filter_info, tumor)
        filter_format.to_csv(f'{os.path.dirname(mut_paths[tumor])}/{tumor}.tmb_mutations.vcf', sep='\t', index=False)
        tmb_num = filter_format.shape[0]
        tmb_value = tmb_num/39 if chip == "wes" else tmb_num/1.9
        result_dic[tumor] = f'{tmb_value:.2f}'
    with open(f'{os.path.abspath(out)}/tmb_result_{chip}.txt', 'w') as result_file:
        result_file.write('sample\ttmb_value\n')
        # 按键排序sorted(result_dic.keys())
        for tumor, tmb_value in sorted(result_dic.items(), key=lambda x: float(x[1]), reverse=True):
            print(f'{tumor} ----- {tmb_value}')
            result_file.write(f'{tumor}\t{tmb_value}\n')


def get_args():
    parser = argparse.ArgumentParser(description='Calculate TMB value by sentieon vcf file.')
    parser.add_argument('--glob_path', '-g', help='注释文件路径，支持"*"', required=True)
    parser.add_argument('--out', '-o', help='输出文件路径', required=True)
    parser.add_argument('--chip', '-c', help='wes按39mb计算，980panel按1.9mb计算', choices=['wes', '980'])
    return parser.parse_args()


if __name__ == "__main__":
    pd.set_option('max_rows', 5)
    pd.set_option('max_columns', 15)
    args = get_args()
    calculate_tmb(args.glob_path, args.chip, args.out)
