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
    print(date[0])
    # vaf_sum = VafSummary(args, date)
    # vaf_sum.sum()


def type_sum(args):
    date = format_date(args.year, args.month, args.day)
    print(date[0])
    # type_sum = TypeSummary(args, date)
    # type_sum.sum()


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
    type_parser.add_argument('--type', '-t', help='变异类型，格式为C>T（包含G>A）', required=True)
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
        self.mutdf = pd.DataFrame()
        self.samples = 0
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
            self.mutdf.to_csv(f'{self.out}/{self.product}_{self.chip}.{self.start_date}_{self.end_date}.somatic_stats.xls', sep='\t', index=False)
        else:
            self.mutdf.to_csv(f'{self.out}/{self.chip}.{self.start_date}_{self.end_date}.somatic_stats.xls', sep='\t', index=False)


class TypeSummary(object):
    def __init__(self, args, date):
        self.chip = args.chip
        self.product = args.product
        self.sumdir = args.sumdir
        self.type = args.type
        if not os.path.exists(args.out): os.mkdir(args.out)
        self.out = args.out
        self.start_date = date[0]
        self.end_date = date[1]
        self.mutdf = pd.DataFrame()
        self.samples = 0

    def sum(self):
        for adir in os.listdir(self.sumdir):
            if not os.path.isdir(f'{self.sumdir}/{adir}'):
                continue
            a_split = adir.split('_')
            if a_split[1].startswith('A') and a_split[3] == self.chip and adir.endswith('auto'):
                if self.product and a_split[2] != self.product:
                    continue
                else:
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
            self.mutdf.to_csv(f'{self.out}/{self.product}_{self.chip}.{self.start_date}_{self.end_date}.somatic_stats.xls', sep='\t', index=False)
        else:
            self.mutdf.to_csv(f'{self.out}/{self.chip}.{self.start_date}_{self.end_date}.somatic_stats.xls', sep='\t', index=False)


if __name__ == "__main__":
    args = get_args()
    args.func(args)
