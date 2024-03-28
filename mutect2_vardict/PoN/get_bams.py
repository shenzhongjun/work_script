#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
根据订单号和样本号获得bam文件路径
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-3-26 09:47:42"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import subprocess


parser = argparse.ArgumentParser(description='根据订单号和样本号获得bam文件路径。')
parser.add_argument('--input', '-i', help='输入文件', required=True)
parser.add_argument('--output', '-o', help='输出文件', required=True)
args = parser.parse_args()
with open(args.input) as f, open(args.output, 'w') as w:
    for line in f:
        chip = line.strip()
        with open(f'{chip}/{chip}samples.txt') as samples:
            for i in samples:
                order = i.strip().split('\t')[5]
                sample = i.strip().split('\t')[2]
                ff_result = subprocess.getoutput(f'ff {order}').strip()
                if ff_result:
                    for order_path in ff_result.split('\n'):
                        ls_result = subprocess.getoutput(f'ls {order_path}/DNA*/bwa_alignment/{sample}/{sample}.final.bam').strip().split('\n')[0]
                        if not ls_result.startswith('ls'):
                            print(f'{chip}\t{sample}\t{order}\t{ls_result}')
                            w.write(f'{chip}\t{sample}\t{order}\t{ls_result}\n')
                            break

