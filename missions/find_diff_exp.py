#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
以一列订单号为输入文件，查找现存目录下MKRN1差异表达信息
"""

import re
import argparse
import subprocess
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='以一列订单号为输入文件，查找现存目录下MKRN1差异表达信息')
    parser.add_argument('--orders', help='订单列表', required=True)
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    argv = get_args()
    with open(argv.orders) as orders, open(argv.out, 'w') as out:
        out.write('订单号\t基因名\t标准化后的序列支持数:检测样本\t标准化后的序列支持数:SRR1691645\t标准化后的序列支持数:SRR1691646\t标准化后的序列支持数:SRR1691647\t标准化后的序列支持数:SRR1691648\tFPKM值:Sample\tFPKM值:Control\tFDR值\tlog2(FoldChange)值\t基因上下调\t是否差异显著基因\n')
        for order in orders:
            order = order.strip()
            order_ff = subprocess.getoutput(f'ff {order}')
            if order_ff:
                for order_path in order_ff.split('\n'):
                    rna_path = subprocess.getoutput(f"find {order_path}/RNA_NBBR/DiffExpression/Sample_vs_Control -type f -name 'Sample_vs_Control.final.txt'")
                    if not rna_path.startswith('find'):
                        with open(rna_path) as exp:
                            for line in exp:
                                if line.startswith('MKRN1'):
                                    out.write(f'{order}\t{line.strip()}\n')
                                    break
                        break

