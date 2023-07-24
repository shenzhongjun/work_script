#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
特定芯片历史订单中变异情况统计。
功能：
1. 统计给定日期区间内某一变异在历史样本中的出现频率（出现样本数/总样本数）
2. 统计给定日期区间内某一变异种类在历史样本中的出现比例（特定变异出现次数/总变异数）
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-4-27 14:30:57"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import sys
from datetime import datetime

import pandas as pd


def vaf_sum(args):
    date = format_date(args.year, args.month, args.day)
    vaf_sum = VafSummary(args, date)
    vaf_sum.sum()


def type_sum(args):
    date = format_date(args.year, args.month, args.day)
    type_sum = TypeSummary(args, date)
    type_sum.sum()


def get_args():
    parser = argparse.ArgumentParser(description='特定芯片或产品历史订单中变异情况统计。')
    subparsers = parser.add_subparsers(title='可选功能')  # , description='统计aaa', help='统计AAA'

    vaf_parser = subparsers.add_parser(name='vafsum', help='统计变异在历史样子中出现频率')
    vaf_parser.add_argument('--chip', '-c', help='芯片类型', required=True)
    vaf_parser.add_argument('--product', '-p', help='产品类型')
    vaf_parser.add_argument('--sumdir', '-s', help='所有样本结果汇总目录',
                            default='/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/')
    vaf_parser.add_argument('--out', '-o', help='芯片变异统计结果输出路径',
                            default=f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/log')
    vaf_parser.add_argument('--day', '-d', help='要统计的日期，支持1-31格式', default='1-31')
    vaf_parser.add_argument('--month', '-m', help='要统计的月份，支持1-3格式', default='1-12')
    vaf_parser.add_argument('--year', '-y', help='要统计的年份，支持2020-2022格式', default='2018-2118')
    vaf_parser.set_defaults(func=vaf_sum)

    type_parser = subparsers.add_parser(name='typesum', help='统计历史特定碱基变异形式比例')
    type_parser.add_argument('--type', '-t', help='变异类型，格式为"C>T"（包含G>A）', required=True)
    type_parser.add_argument('--chip', '-c', help='芯片类型', required=True)
    type_parser.add_argument('--product', '-p', help='产品类型')
    type_parser.add_argument('--sumdir', '-s', help='所有样本结果汇总目录',
                             default='/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/')
    type_parser.add_argument('--out', '-o', help='芯片变异统计结果输出路径',
                             default=f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/log')
    type_parser.add_argument('--day', '-d', help='要统计的日期，支持1-31格式', default='1-31')
    type_parser.add_argument('--month', '-m', help='要统计的月份，支持1-3格式', default='1-12')
    type_parser.add_argument('--year', '-y', help='要统计的年份，支持2020-2022格式', default='2018-2118')
    type_parser.set_defaults(func=type_sum)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        return parser.parse_args()


def format_date(year, month, day):
    """处理参数输入的日期数字，返回开始日期和截止日期"""

    def check_date(date, date_type):
        date = int(date)
        if date_type == 'd':
            assert 1 <= date <= 31, '请检查日期输入！'
        elif date_type == 'm':
            assert 1 <= date <= 12, '请检查日期输入！'
        return date

    def get_date_list(date, date_type):
        date_list = []
        if '-' in date:
            start, end = date.split('-')
            for l in range(int(start), int(end) + 1):
                date_list.append(check_date(l, date_type))
        else:
            date_list.append(check_date(date, date_type))
        return date_list

    year_list = sorted(get_date_list(year, 'y'))
    month_list = sorted(get_date_list(month, 'm'))
    day_list = sorted(get_date_list(day, 'd'))
    start_date = datetime(year_list[0], month_list[0], day_list[0])
    # 防止日期上限（31）大于月份实际上限（2月）
    while True:
        try:
            end_date = datetime(year_list[-1], month_list[-1], day_list[-1])
        except ValueError:
            day_list.pop()
        else:
            break
    return start_date, end_date


class Summary(object):
    def __init__(self, args, date):
        self.chip = args.chip
        self.product = args.product
        self.sumdir = args.sumdir
        if not os.path.exists(args.out): os.mkdir(args.out)
        self.out = args.out
        self.start_date = date[0]
        self.end_date = date[1]
        self.dirs = self.get_dirs()

    def get_dirs(self):
        dirs = []
        for adir in os.listdir(self.sumdir):
            if not os.path.isdir(f'{self.sumdir}/{adir}'):
                continue
            a_split = adir.split('_')
            if a_split[1].startswith('A') and a_split[3] == self.chip and adir.endswith('auto'):
                if self.product and a_split[2] != self.product:
                    continue
                else:
                    dirs.append(f'{self.sumdir}/{adir}')
        return dirs


class VafSummary(Summary):
    def __init__(self, args, date):
        super().__init__(args, date)
        self.samples = 0
        self.mutdf = pd.DataFrame()

    def sum(self):
        for adir in self.dirs:
            for i in os.listdir(f'{self.sumdir}/{adir}'):
                merge_file = f'{self.sumdir}/{adir}/{i}'
                ctime = datetime.fromtimestamp(os.stat(merge_file).st_ctime)
                if i.endswith('_Somatic.xlsx') and self.start_date < ctime < self.end_date:
                    df = pd.read_excel(merge_file, usecols='A, H, J, L, M, R:T, W, X, Z:AH, CC')
                    df = df.drop_duplicates(subset='突变编号', keep='first', ignore_index=True)
                    df = df[(~df['结论'].isin(['遗传性突变'])) &
                            (df['突变可靠性\n'] == '相对可靠') &
                            (~df['突变类型'].isin(['同义突变', '非编码区突变', '剪切位点附近的突变']))]
                    self.mutdf = self.mutdf.append(df, ignore_index=True)
                    self.samples += 1
        mut_count_df = pd.DataFrame({'突变编号': self.mutdf['突变编号'].value_counts().index,
                                     '该突变检出的个数': self.mutdf['突变编号'].value_counts()})
        mut_count_df = mut_count_df[mut_count_df['该突变检出的个数'] >= 3]  # 只统计样本数大于2个的变异
        mut_join_df = self.mutdf.groupby('突变编号')['结论'].apply(lambda x: ';'.join(x.unique())).to_frame()  # 不同结论合并
        self.mutdf = self.mutdf.drop_duplicates(subset='突变编号', keep='first', ignore_index=True).drop('结论', axis=1)
        self.mutdf = self.mutdf.merge(mut_count_df, on='突变编号')
        self.mutdf = self.mutdf.merge(mut_join_df, on='突变编号')
        self.mutdf['样本数'] = self.samples
        self.mutdf['频率'] = self.mutdf['该突变检出的个数'] / self.mutdf['样本数']
        self.mutdf = self.mutdf.reindex(['突变编号', '结论', '基因', 'mRNA变体名称(NM)', '染色体名称', '核酸改变',
                                         '氨基酸改变', '突变类型', '生物学意义', 'rs号', '次要基因型频率',
                                         '1000genomeMAF', 'CHBmaf', 'CHSmaf', 'ExAC-eas', 'ExAC-sas', 'ExAC-biggest',
                                         'meMAF_nom', 'ESP_MAF', '该突变检出的个数', '样本数', '频率'], axis=1)
        self.mutdf = self.mutdf.round({'频率': 4}).sort_values(['突变编号'])
        if self.product:
            self.mutdf.to_csv(f'{self.out}/{self.product}_{self.chip}.{self.start_date.date()}_{self.end_date.date()}.somatic_stats.xls',
                              sep='\t', index=False)
        else:
            self.mutdf.to_csv(f'{self.out}/{self.chip}.{self.start_date.date()}_{self.end_date.date()}.somatic_stats.xls',
                              sep='\t', index=False)


class TypeSummary(Summary):
    def __init__(self, args, date):
        super().__init__(args, date)
        self.mut_type = args.type.split('>')
        trantab = args.type.maketrans('ATCG', 'TAGC')
        self.mut_type_flip = args.type.translate(trantab).split('>')
        print(f'\n要统计的突变类型为[{">".join(self.mut_type)}]及[{">".join(self.mut_type_flip)}]\n')
        self.mut_num = 0
        self.mut_type_num = 0
        
    def sum(self):
        for adir in self.dirs:
            for i in os.listdir(adir):
                merge_file = f'{adir}/{i}'
                ctime = datetime.fromtimestamp(os.stat(merge_file).st_ctime)
                if i.endswith('_Somatic.xlsx') and self.start_date < ctime < self.end_date:
                    df = pd.read_excel(merge_file, usecols='A, P, Q')
                    df = df.drop_duplicates(subset='突变编号', keep='first', ignore_index=True)
                    self.mut_num += df.shape[0]
                    df = df[((df['参考碱基'] == self.mut_type[0]) & (df['突变碱基'] == self.mut_type[1])) |
                            ((df['参考碱基'] == self.mut_type_flip[0]) & (df['突变碱基'] == self.mut_type_flip[1]))]

                    self.mut_type_num += df.shape[0]

        if self.product:
            writer = open(f'{self.out}/{self.product}_{self.chip}.{self.start_date.date()}_{self.end_date.date()}.{"-".join(self.mut_type)}_stats.xls', 'w')
        else:
            writer = open(f'{self.out}/{self.chip}.{self.start_date.date()}_{self.end_date.date()}.{"-".join(self.mut_type)}_stats.xls', 'w')

        print(f'给定时间段内[{self.chip}]芯片变异共：{self.mut_num}个\n')
        print(f'统计{">".join(self.mut_type)}及{">".join(self.mut_type_flip)}类型突变共：{self.mut_type_num}个，'
              f'占比为{self.mut_type_num/self.mut_num*100:.2f}%')

        writer.write('芯片类型\t突变类型\t突变数目\t总突变数目\t突变比例(%)\n')
        writer.write(f'{self.chip}\t{">".join(self.mut_type)}({">".join(self.mut_type_flip)})\t{self.mut_type_num}\t'
                     f'{self.mut_num}\t{self.mut_type_num/self.mut_num*100:.2f}\n')


if __name__ == "__main__":
    args = get_args()
    args.func(args)
