#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
对给定样本进行质控统计。
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-6-29 10:21:44"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import re
import sys
import glob
from datetime import datetime
import time
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='对给定样本进行质控统计。')
    parser.add_argument('--sumdir', '-s', help='所有样本结果汇总目录，默认为/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/',
                        default='/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/')
    parser.add_argument('--out', '-o', help='样本质控结果输出路径，默认为脚本目录下log文件夹',
                        default=f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/log')
    parser.add_argument('--list', '-l', help='要统计的样本列表', default='sample_list')
    return parser.parse_args()


class QcSummary(object):
    def __init__(self, args, sample_list, batch_list):
        self.sample_list = sample_list
        self.batch_list = batch_list
        self.sumdir = args.sumdir
        if not os.path.exists(args.out): os.mkdir(args.out)
        self.out = args.out
        self.sampledf = pd.DataFrame()

        for batch, sample in zip(self.batch_list, self.sample_list):
            glob_path = glob.glob(fr'{self.sumdir}/*{sample}*{batch}*/{sample}_*.xlsx')
            if not glob_path:
                sample = sample.strip('DZ')
                glob_path = glob.glob(fr'{self.sumdir}/*{sample}*{batch}*/{sample}_*.xlsx')
                if not glob_path:
                    print(f'不存在{sample}文件夹或{sample}文件夹下没有质控表！')
                    continue
            elif len(glob_path) > 1:
                print(f'存在多个{sample}文件夹或{sample}文件夹下有多个质控表！使用第一个进行统计！')
                print(glob_path)

            qc_file = glob_path[0]
            try:
                df = pd.read_excel(qc_file)
            except Exception as ex:
                print(f'{sample} qc文件非Excel格式，尝试使用txt格式打开.{ex}')
                df = pd.read_table(qc_file)  # 有可能是0.9流程生成的文本文件但是后缀为xlsx

            self.sampledf = self.sampledf.append(df, ignore_index=True)
        self.sampledf = self.sampledf.drop_duplicates(subset='样本编号', keep='first', ignore_index=True)
        self.sampledf.to_csv(f'{self.out}/{time.strftime("%Y%m%d")}.samples.qc_merge.txt', sep='\t', index=False)


if __name__ == "__main__":
    args = get_args()
    list_path = args.list
    sample_list = []
    batch_list = []
    with open(list_path) as f:
        for i in f:
            batch_list.append(i.strip().split('\t')[0])
            sample_list.append(i.strip().split('\t')[1])
    qc_sum = QcSummary(args, sample_list, batch_list)
