#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
定期批量样本质控统计。
功能：
1. 识别生产样本：样本以A开头，文件夹以auto结尾
2. 用ctime判断文件修改时间，符合所给时间区间的样本合并质控表，以样本名去重
3. 根据规则进行质控，结果输出到文件
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-4-2 09:02:57"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import re
import sys
from datetime import datetime

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='定期批量样本质控统计。')
    parser.add_argument('--sumdir', '-s', help='所有样本结果汇总目录，默认为/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/',
                        default='/mnt/share03/tumor_pipeline/Somatic/DNA/Other/New_report/')
    parser.add_argument('--out', '-o', help='样本质控结果输出路径，默认为脚本目录下log文件夹',
                        default=f'{os.path.split(os.path.realpath(sys.argv[0]))[0]}/log')
    parser.add_argument('--day', '-d', help='要统计的日期，支持1-31格式，默认为1-31', default='1-31')
    parser.add_argument('--month', '-m', help='要统计的月份，支持1-3格式，必需参数', required=True)
    parser.add_argument('--year', '-y', help='要统计的年份，支持2020-2022格式，默认为2023', default='2023')
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


class QcSummary(object):
    def __init__(self, args, date):
        self.chip_type = {
            'large_panel': ['NPC980', 'NST', 'NLP'],
            'small_panel': ['NPC116', 'NPC69'],
            'wes': ['NT01T', 'IDT'],
            'rna_panel': ['NSRF', 'NBRF'],
            'rna_seq': ['mRNA'],
            'immune': ['PCR'],
            'other': ['NVT', 'NCDT', 'NMM26', 'NHLA', 'NCP1']
        }
        self.sumdir = args.sumdir
        if not os.path.exists(args.out): os.mkdir(args.out)
        self.out = args.out
        self.start_date = date[0]
        self.end_date = date[1]
        self.sampledf = pd.DataFrame()
        self.samples = set()  # 一个set用来装本期所有样本对象

        for adir in os.listdir(self.sumdir):
            if not os.path.isdir(f'{self.sumdir}/{adir}'):
                continue
            if adir.split('_')[1].startswith('A') and adir.endswith('auto'):    # 生产订单特征，与研发项目区分
                for i in os.listdir(f'{self.sumdir}/{adir}'):
                    qc_file = f'{self.sumdir}/{adir}/{i}'
                    ctime = datetime.fromtimestamp(os.stat(qc_file).st_ctime)
                    if i.endswith('rmdup_QC.xlsx') and self.start_date < ctime < self.end_date:
                        try:
                            df = pd.read_excel(qc_file)
                        except:
                            df = pd.read_table(qc_file)     # 有可能是0.9流程生成的文本文件但是后缀为xlsx
                        self.sampledf = self.sampledf.append(df, ignore_index=True)
                        for l in range(0, 2):
                            sample = Sample(df.iloc[l])
                            sample.order = re.search('(DDN.*?)_', adir).group(1)
                            sample.chip_type = self.chip_type
                            self.samples.add(sample)
        self.sampledf = self.sampledf.drop_duplicates(subset='样本编号', keep='first', ignore_index=True)
        self.sampledf.to_csv(f'{self.out}/{self.start_date.date()}_to_{self.end_date.date()}.qc_merge.txt',
                             sep='\t', index=False)
        self.sampledf.to_excel(f'{self.out}/{self.start_date.date()}_to_{self.end_date.date()}.qc_merge.xlsx',
                               index=False)
        # print(self.sampledf)

    def qc_summary(self):
        sum_dict = {
            'WES': {'all_sample': 0, 'nopass': 0},
            '大panel': {'all_sample': 0, 'nopass': 0},
            '小panel': {'all_sample': 0, 'nopass': 0},
            '其它': {'all_sample': 0, 'nopass': 0},
            '汇总': {'all_sample': 0, 'nopass': 0}
        }
        sum_list = []   # 不合格统计
        nopass_list = []    # 不合格明细
        sum_writer = pd.ExcelWriter(
            f'{self.out}/{self.start_date.date()}_to_{self.end_date.date()}.下机数据质控失败率统计.xlsx')

        for sample in self.samples:
            sample_type = sample.type.replace('wes', 'WES').replace('large_panel', '大panel'). \
                replace('small_panel', '小panel').replace('other', '其它')
            sum_dict[sample_type]['all_sample'] += 1
            sum_dict['汇总']['all_sample'] += 1
            if not sample.passed:
                sum_dict[sample_type]['nopass'] += 1
                sum_dict['汇总']['nopass'] += 1
                nopass_list.append([sample.order, sample.name, sample.chip, sample.unpass_reason, sample.unpass_info])

        for i in sum_dict:
            if sum_dict[i]["all_sample"] > 0 and sum_dict[i]["nopass"] > 0:
                nopass_rate = f'{(sum_dict[i]["nopass"] / sum_dict[i]["all_sample"]) * 100:.2f}%'
            else:
                nopass_rate = '0'
            sum_list.append([i, sum_dict[i]["all_sample"], sum_dict[i]["nopass"], nopass_rate])

        sum_df = pd.DataFrame.from_records(
            sum_list, columns=['产品类型', '下机样本总数', '质控不合格数', '不合格率'])
        sum_df.to_excel(sum_writer, sheet_name="质控失败率统计", index=False)
        print(f'质控失败率统计：\n{sum_df}')

        sum_df2 = pd.DataFrame.from_records(
            nopass_list, columns=['订单编号', '样本号', '芯片类型', '不合格原因', '不合格明细'])
        sum_df2.to_excel(sum_writer, sheet_name="质控不合格明细", index=False)
        sum_writer.save()
        print(f'质控不合格明细：\n{sum_df2}')


class Sample(object):
    def __init__(self, df):
        self.chip_type = {}
        self.order = ''
        self.name = df['样本编号']
        self.batch = df['批次号']
        self.product = df['项目代码']
        self.chip = df['捕获芯片']
        self.captive_rate = df['捕获率']
        self.coverage = df['目标区域覆盖度(%)']
        self.data_size = df['测序数据量(M)']
        self.dup = df['重复率']
        self.deepth = df['去重复深度']
        self.insert_size = df['插入片段均值']
        self.q30 = df['Q30']
        # self.coverage_300x = df['<=300x目标区域覆盖度(%)']
        # self.coverage_600x = df['<=600x目标区域覆盖度(%)']
        # self.alignment_rate = df['参考序列比对率']
        # self.paired_consistency = df['配对样本一致性']
        self.is_normal = True if 'DZ' in self.name else False
        self.unpass_dict = {'reason': [], 'info': []}
        self.unpass_reason = ''
        self.unpass_info = ''

    @property
    def type(self):
        return [i for i in self.chip_type if self.chip in self.chip_type[i]][0]

    @property
    def passed(self):
        if self.filter_deepth() and self.filter_q30() and self.filter_insert_size() and self.filter_dup():
            return True
        else:
            self.unpass_reason = ';'.join(self.unpass_dict['reason'])
            self.unpass_info = ';'.join(self.unpass_dict['info'])
            return False

    def filter_q30(self):
        if float(self.q30) <= 85:
            self.unpass_dict['reason'].append('Q30低')
            self.unpass_dict['info'].append(str(self.q30))
            return False
        else:
            return True

    def filter_dup(self):
        if float(self.dup) >= 0.9:
            self.unpass_dict['reason'].append('重复率异常')
            self.unpass_dict['info'].append(str(self.dup))
            return False
        else:
            return True

    def filter_deepth(self):  # RNA和免疫组库未实装
        if self.is_normal:
            if (int(self.deepth) <= 150 and self.type in ['large_panel', 'small_panel', 'other'] or
                    int(self.deepth) <= 70 and self.type == 'wes'):
                self.unpass_dict['reason'].append('对照有效深度不足')
                self.unpass_dict['info'].append(str(self.deepth))
                return False
            else:
                return True
        else:
            solid_wes = (self.product == 'Ncet' or self.product == 'Ncec') and self.deepth <= 300 and self.type == 'wes'
            other_wes = (self.product != 'Ncet' and self.product != 'Ncec') and self.deepth <= 450 and self.type == 'wes'
            if (solid_wes or other_wes or
                    (int(self.deepth) <= 800 and (self.type in ['large_panel', 'small_panel'] or
                                                  (self.type == 'other' and self.chip in ['NVT', 'NCDT', 'NMM26', 'NCP1']))) or
                    (int(self.deepth) <= 200 and self.chip == 'NHLA')):
                self.unpass_dict['reason'].append('肿瘤有效深度不足')
                self.unpass_dict['info'].append(str(self.deepth))
                return False
            else:
                return True

    def filter_insert_size(self):
        if ((self.type in ['wes', 'large_panel', 'small_paned'] or self.chip in ['NVT', 'NCDT', 'NMM26', 'NCP1'])
                and int(self.insert_size) <= 80):
            self.unpass_dict['reason'].append('插入片段长度不足')
            self.unpass_dict['info'].append(str(self.insert_size))
            return False
        else:
            return True

    def filter_data_size(self):  # RNA和免疫组库未实装，DNA仅有风险指标，无不合格指标
        pass

    def filter_alignment_rate(self):  # QC表无此项目，未实装。均一性质控、配对样本一致性也没有。
        pass


if __name__ == "__main__":
    args = get_args()
    date = format_date(args.year, args.month, args.day)
    qc_sum = QcSummary(args, date)
    qc_sum.qc_summary()
