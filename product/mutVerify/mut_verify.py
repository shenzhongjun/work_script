#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对ddPCR引物探针目录.xlsx进行数据清洗
对输入的位点是否能进行ddPCR实验验证进行判断
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-7-20 17:17:19"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import sys
import re
import pandas as pd


def get_mut(x):
    if 'MET Splice' in x:
        gene = 'MET'
        change = re.search(r'Mutation\((c\.[0-9]+.*)\)', x).group(1)
        change = change.strip('c.')
    else:
        gene, change = x.split(' ')[1:3]
        change = change.strip('p.').strip('c.')
    return [gene, change]


if __name__ == "__main__":
    argv = sys.argv[1]
    pd.set_option('display.max_colwidth', 100)
    if argv.endswith('xlsx'):
        df = pd.read_excel(argv)
        df = df.fillna(method='ffill', axis=0)
        df = df[~df['Mutation'].str.contains('WT')]
        df = df[df['Application Type:'] == 'Mu-DdPCR Package']
        df[['gene', 'change']] = df['Mutation'].apply(get_mut).apply(pd.Series)
        df['突变'] = df['gene'].str.cat(df['change'], sep=':')
        df = df.rename({'配套出售货号:': '验证试剂盒货号'}, axis=1)
        df = df[['Mutation', '突变', '验证试剂盒货号']]
        df.to_csv('ddPCR.clean.txt', sep='\t', index=False)
        print(df)
    else:
        gene, change = argv.split(':')
        change = change.strip('p.').strip('c.')
        argv = ':'.join([gene, change])
        df = pd.read_table(f'{os.path.dirname(os.path.realpath(__file__))}/ddPCR.clean.txt')
        result = df[df['突变'] == argv]
        if result.shape[0] > 0:
            print(f"该位点可以进行实验验证：\n{result.drop('突变', axis=1).to_string(index=False)}")
        else:
            print('该位点无法进行实验验证')

