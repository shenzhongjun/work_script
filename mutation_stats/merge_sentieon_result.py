#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Sentieon call变异结果添加到merge表最后；
2022年1月24日：经讨论决定不添加
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2021-12-17 16:33:39"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import io
import gzip
import argparse
import pandas as pd
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description='Sentieon call变异结果添加到merge表最后')
    parser.add_argument('--merge_table', '-m', help='txt格式merge表路径', required=True)
    parser.add_argument('--sentieon_vcf', '-v', help='sentieon *.PASS.vcf 路径', required=True)
    parser.add_argument('--tumor_sample', '-t', help='肿瘤样本名', required=True)
    parser.add_argument('--out_prefix', '-o', help='输出文件前缀', required=True)
    return parser.parse_args()


def trans_pos(x):
    """
    变异类型 |    标准vcf变异格式     | 老流程变异格式
    del      | chr6_157099871_TCCG_T  | chr6_157099872_CCG_-
    insert   | chr6_31324205_G_GCC to | chr6_31324205_-_CC
    因使用chr_position合并，只处理del的坐标即可
    """
    chrom, pos, ref = x.split('_')[0:3]
    pos = int(pos)
    if len(ref) > 1:
        pos += 1
    return '_'.join([chrom, str(pos)])


def get_sentieon_df(file_path, tumor_sample):
    """
    related header lines:
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">

    vcf FORMAT：
    GT:AD:AF:AFDP:ALTHC:ALT_F1R2:......

    vcf TUMOR：
    0/1:26151,83:0.003:26234:83:83:......
    """
    tmp = io.StringIO()
    if 'gz' in file_path:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp.write(line)
    else:
        with open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp.write(line)
    df = pd.read_table(io.StringIO(tmp.getvalue()), dtype={'POS': 'str'})   # 注：加comment='#'参数即可实现跳过'##'行！

    df['Sentieon_突变位置'] = df['#CHROM'].str.cat([df['POS'], df['REF'], df['ALT']], sep='_')
    df['Position'] = df['Sentieon_突变位置'].apply(trans_pos)  # 处理成老流程坐标且只保留chr_position
    df['Sentieon_参考深度'] = df[tumor_sample].str.split(':', expand=True)[1].str.split(',', expand=True)[0]
    df['Sentieon_突变深度'] = df[tumor_sample].str.split(':', expand=True)[1].str.split(',', expand=True)[1]
    df['Sentieon_总深度'] = df['Sentieon_参考深度'].astype(np.int64) + df['Sentieon_突变深度'].astype(np.int64)
    df['Sentieon_突变率'] = df[tumor_sample].str.split(':', expand=True)[2]
    df['Sentieon_总深度'] = df['Sentieon_总深度'].astype('str')
    df = df[['Position', 'Sentieon_突变位置', 'Sentieon_总深度', 'Sentieon_突变深度', 'Sentieon_突变率']]
    return df


if __name__ == "__main__":
    args = get_args()
    merge_table = pd.read_table(os.path.abspath(args.merge_table))
    merge_table['Position'] = merge_table['染色体名称'].str.cat(merge_table['起始位置'].astype('str'), sep='_')
    sentieon_df = get_sentieon_df(os.path.abspath(args.sentieon_vcf), args.tumor_sample)
    merged_df = merge_table.merge(sentieon_df, how='left', on='Position')
    merged_df_sen_uniq = merge_table.merge(sentieon_df, how='right', on='Position')
    merged_df_sen_uniq = merged_df_sen_uniq[merged_df_sen_uniq['突变编号'].isna()]
    merged_df_sen_uniq = merged_df_sen_uniq[['Sentieon_突变位置', 'Sentieon_总深度', 'Sentieon_突变深度', 'Sentieon_突变率']]
    merged_df.drop(columns=['Position', 'Sentieon_突变位置']).to_csv(f'{args.out_prefix}.txt', sep='\t', index=False)
    merged_df.drop(columns=['Position', 'Sentieon_突变位置']).to_excel(f'{args.out_prefix}.xlsx', index=False)
    # merged_df_sen_uniq.to_csv(f'{args.out_prefix}.sentieon_uniq.txt', sep='\t', index=False)
