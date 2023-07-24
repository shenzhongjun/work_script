#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Calculate TMB value by annotated vcf file.
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-9-7 10:24:41 "
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import argparse
import pandas as pd
import numpy as np


def filter_by_info(df):
    """
    根据vcf文件INFO列过滤
    """
    df = df[df['INFO'].str.contains('SOMATIC;SS=2')]    # filter SOMATIC
    df = df[df['INFO'].str.split(';', expand=True)[3].str.split('=', expand=True)[1].astype(np.int64) >= 200]   # filter SOMATIC by SSC
    filtered_df = df[df['INFO'].str.contains('nonsynonymous|stopgain|stoploss|nonframeshift')]
    # print(filtered_df)
    return filtered_df


def filter_by_format(df):
    """
    根据vcf文件FORMAT列过滤
    """
    df_n = df['NORMAL'].str.split(':', expand=True)
    df_t = df['TUMOR'].str.split(':', expand=True)
    df_n.columns = ['GT', 'GQ', 'DP', 'RD', 'AD', 'AF', 'DP4']
    df_t.columns = ['GT', 'GQ', 'DP', 'RD', 'AD', 'AF', 'DP4']
    # df_n = df_n[(df_n['DP'].astype(np.int64) > 50) &
    #             (df_n['DP4'].str.split(',', expand=True)[2].astype(np.int64) +
    #              df_n['DP4'].str.split(',', expand=True)[3].astype(np.int64) >= 7)]
    # df_n = df_n[(df_n['DP'].astype(np.int64) > 50) &
    #             (df_n['AF'].str.split('%', expand=True)[0].astype(np.float64) <= 0.35)]
    df_n = df_n[(df_n['DP'].astype(np.int64) > 50)]
    print(df_n)
    df_t = df_t[(df_t['DP'].astype(np.int64) > 50) &
                (df_t['AF'].str.split('%', expand=True)[0].astype(np.float64) >= 5) &
                (df_t['DP4'].str.split(',', expand=True)[2].astype(np.int64) +
                 df_t['DP4'].str.split(',', expand=True)[3].astype(np.int64) >= 8)]
    print(df_t)
    df_n_t = pd.concat([df_n, df_t], axis=1, join='inner')
    df_final = df.loc[df_n_t.index.tolist()]
    return df_final


def calculate_tmb(obj):
    raw_df = pd.read_table(obj, header=30)
    filter_info = filter_by_info(raw_df)
    filter_format = filter_by_format(filter_info)
    # print(filter_info, filter_format, sep='\n')
    path = os.path.split(os.path.realpath(args.anno))[0]
    tmb_name = f'{os.path.splitext(os.path.basename(args.anno))[0]}.SSC_200.TMB.vcf'
    filter_format.to_csv(os.path.join(path, tmb_name), sep='\t', index=False)

    tmb_num = filter_format.shape[0]
    TMB_value = f'{tmb_num / 39:.2f}'
    stat = pd.DataFrame({'TMB_somatic_count': [tmb_num], 'TMB': [TMB_value]})
    stat_name = f'{os.path.splitext(os.path.basename(args.anno))[0]}.SSC_200.TMB_stat.xls'
    stat.to_csv(os.path.join(path, stat_name), sep='\t', index=False)
    print(stat)


def get_args():
    parser = argparse.ArgumentParser(description='Calculate TMB value by annotated vcf file.')
    parser.add_argument('--anno', help='the annotations file', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    pd.set_option('max_rows', 5)
    pd.set_option('max_columns', 15)
    args = get_args()
    calculate_tmb(args.anno)
