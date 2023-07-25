#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Annovar txt注释结果过滤基因列表，用于一项科研任务
"""

import argparse
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='bed与vcf匹配，如果在vcf中，修改bed并加一列check，否则加nocheck')
    parser.add_argument('--genes', help='基因列表文件路径', required=True)
    parser.add_argument('--txt', help='Annovar txt注释文件路径', required=True)
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    pd.set_option('max_columns', 15)
    argv = get_args()
    genes = []
    with open(argv.genes) as f:
        for line in f:
            genes.append(line.strip())
    df = pd.read_table(argv.txt)
    df = df[df['Gene.refGene'].isin(genes)]
    df.to_csv(argv.out, sep='\t', index=False)

