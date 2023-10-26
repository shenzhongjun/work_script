#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Annovar txt注释结果过滤基因列表，用于一项科研任务
"""

import re
import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='bed与vcf匹配，如果在vcf中，修改bed并加一列check，否则加nocheck')
    parser.add_argument('--genes', help='基因列表文件路径', required=True)
    parser.add_argument('--txt', help='Annovar txt注释文件路径', required=True)
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


def get_vaf(x):
    dp = re.search('DP=(.*?);', x).group(1)
    vd = re.search('VD=(.*?);', x).group(1)
    if int(vd) > int(dp):
        vd = dp
    vaf = f'{int(vd)/int(dp):.3f}'
    return dp, vd, vaf


if __name__ == "__main__":
    pd.set_option('max_columns', 15)
    argv = get_args()
    genes = []
    with open(argv.genes) as f:
        for line in f:
            genes.append(line.strip())
    df = pd.read_table(argv.txt)
    df = df[(df['Gene.refGene'].isin(genes)) & (df['Otherinfo10'] == 'PASS')]
    df[['Otherinfo1', 'Otherinfo2', 'Otherinfo3']] = df['Otherinfo11'].apply(get_vaf).apply(pd.Series)
    df = df.rename({'Otherinfo1': '总深度', 'Otherinfo2': '突变深度', 'Otherinfo3': '突变频率'}, axis=1)
    df.to_csv(argv.out, sep='\t', index=False)

