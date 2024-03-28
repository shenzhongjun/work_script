#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
免疫组库订单获得均一化后的count结果
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-3-4 17:43:26"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import argparse
import subprocess
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='从gff文件中获得两个融合断点的上下游各500bp位置信息，从fa中提取对应序列')
    parser.add_argument('--fusion1', help='融合位点A，例如"chr2:42522656"', required=True)
    parser.add_argument('--fusion2', help='融合位点B，例如"chr2:29446394"', required=True)
    parser.add_argument('--gene1', help='融合基因A', required=True)
    parser.add_argument('--gene2', help='融合基因B', required=True)
    parser.add_argument('--gtf', help='gtf文件路径', required=True)
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    libdf = pd.read_table('factors.txt')
    print(libdf)
    with open('file_list.txt') as file:
        for path in file:
            path = path.strip()
            filename = os.path.splitext(path)[0]
            library = filename.split('/')[-1].split('.')[0]
            factor = libdf[libdf['library'] == library]['factor'].values[0]
            print(library, factor)
            # with open(path) as f, open(f'{filename}.new.xls', 'w') as w:
            rawdf = pd.read_table(path)
            # print(rawdf)
            # rawdf['new_count'] = (rawdf['count'] * factor).astype(int)
            new_count = (rawdf['count'] * factor).astype(int)  # 计算并转换为整数类型
            rawdf.insert(0, 'new_count', new_count)
            print(rawdf)
            rawdf.to_csv(f'{filename}.new.xls', sep='\t', index=False)
