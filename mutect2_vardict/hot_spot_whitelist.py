#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mutect2和Vardict合并后的Vcf文件过热点突变超级白名单
仅供v1.0生产流程使用
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-1-30"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import re
import pandas as pd
import argparse
import gzip
from mutect2_vardict_vcf_merger import Vcf


def get_args():
    parser = argparse.ArgumentParser(description='Vcf文件过热点突变超级白名单。')
    parser.add_argument('--input_vcf', '-i', help='vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf.gz输出文件路径', required=True)
    parser.add_argument('--hot_db', '-d', help='热点突变数据库',
                        default=f'{os.path.dirname(os.path.realpath(__file__))}/hot_spot_whitelist.xlsx')
    return parser.parse_args()


class HotSpotVcf(Vcf):
    def __init__(self, tag, path, input_vcf, hot_db):
        super().__init__(tag, path)
        self.db = hot_db
        self.file_format = input_vcf.file_format
        self.filter = input_vcf.filter
        self.filter.append(
            '##FILTER=<ID=not_hotspot,Description="Detected by Vardict but not in the hot spot whitelist">')
        self.info = input_vcf.info
        self.info.append('##INFO=<ID=HotSpot,Number=0,Type=Flag,Description="Variant is in hot spot whitelist">')
        self.info.append(
            '##INFO=<ID=SuperHotSpot,Number=0,Type=Flag,Description="Variant is in super hot spot whitelist">')
        self.format = input_vcf.format
        self.source = input_vcf.source
        other = f'##Filter_HotSpotdb_Command="python {os.path.abspath(__file__)} --input_vcf {input_vcf.path} --out_vcf {self.path} --hot_db {self.db}'
        self.other = input_vcf.other
        self.other.append(other)
        self.headerline = input_vcf.headerline
        self.contig = input_vcf.contig
        self.mut_df = input_vcf.mut_df
        self.tumor = input_vcf.tumor

    def get_hotspot(self):
        hot_df = pd.read_excel(self.db, names=['基因', '转录本', '位置', '变异', '变异类型', '是否为超级白名单', '备注', '临床意义'])
        hot_df = hot_df[hot_df['变异类型'] != '特殊检测']
        hot_df['AAChange'] = hot_df['基因'].str.cat(hot_df['变异'], sep=':')
        hot_df = hot_df[['AAChange', '是否为超级白名单']]
        hot_df = hot_df.fillna('否')
        hot_dic = dict(zip(hot_df['AAChange'], hot_df['是否为超级白名单']))

        def filterdb(x):
            filter_tags = set(x['filter'].split(';'))
            infos = x['info'].split(';')
            aa_change = infos[-2].split('=')[1]
            if x['format'].startswith('GT:AD:AF:DP'):    # Mutect2格式FORMAT列
                vaf = float(x[self.tumor].split(':')[2])
                reads = int(x[self.tumor].split(':')[1].split(',')[1])
            else:
                vaf = float(x[self.tumor].split(':')[6])
                reads = int(x[self.tumor].split(':')[2])

            if aa_change not in ['.', 'UNKNOWN'] and 'p.' in aa_change:
                # i.split(":")[0]是基因名，i.split(":")[4][2:]是蛋白质改变（去掉前面的p.与白名单进行匹配）
                mut_set = {f'{i.split(":")[0]}:{i.split(":")[4][2:]}' for i in aa_change.split(',') if 'p.' in i}
                hot_mut = mut_set.intersection(set(hot_dic.keys()))
                if hot_mut:
                    if hot_dic[list(hot_mut)[0]] == '是':    # 超级白名单提高报出优先级并加入基础过滤
                        x['info'] = f"{x['info']};SuperHotSpot"
                        if vaf >= 0.004 and reads >= 4 and len(x['ref']) < 50 and len(x['alt']) < 50:
                            filter_tags.add('PASS')
                            if 'germline_risk' in filter_tags:
                                filter_tags.remove('germline_risk')
                    else:   # 普通白名单来源复杂、可靠性较低，视作普通变异，仅添加tag，不再特殊对待
                        x['info'] = f"{x['info']};HotSpot"

            if 'SOFT=Vardict;' in x['info'] and 'PASS' in filter_tags and 'SuperHotSpot' not in x['info'] and 'germline' not in filter_tags:
                filter_tags.add('not_hotspot')
                filter_tags.remove('PASS')

            x['filter'] = ';'.join(sorted(list(filter_tags)))
            return x

        self.mut_df = self.mut_df.apply(filterdb, axis=1)


if __name__ == "__main__":
    args = get_args()
    input_vcf = Vcf('Merged', args.input_vcf)
    input_vcf.get_header()
    input_vcf.get_mutdf()
    output_vcf = HotSpotVcf('HotSpot', args.output_vcf, input_vcf, args.hot_db)
    output_vcf.get_hotspot()
    output_vcf.write_vcf()
