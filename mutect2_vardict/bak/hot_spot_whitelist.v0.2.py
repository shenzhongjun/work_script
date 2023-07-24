#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Vardict Vcf文件筛选级白名单热点突变，兼容特殊变异检测
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-2-6"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import re
import pandas as pd

from handle_vcf import Vcf


def get_args():
    parser = argparse.ArgumentParser(description='Vardict Vcf文件过热点突变超级白名单，兼容特殊变异检测。')
    parser.add_argument('--input_vcf', '-i', help='vcf输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf输出文件路径', required=True)
    parser.add_argument('--hot_db', '-d', help='热点突变数据库',
                        default=f'{os.path.dirname(os.path.realpath(__file__))}/hot_spot_whitelist.xlsx')
    return parser.parse_args()


class HotSpotVcf(Vcf):
    def __init__(self, tag, path, input_vcf, hot_db):
        super().__init__(tag, path)
        self.db = hot_db
        self.file_format = input_vcf.file_format
        self.filter = input_vcf.filter
        self.info = input_vcf.info
        self.info.append(
            '##INFO=<ID=SuperHotSpot,Number=0,Type=Flag,Description="Variant is in super hot spot whitelist.">')
        self.info.append(
            '##INFO=<ID=SuperHotSpot_Special,Number=0,Type=Flag,Description="Variant is in super hot spot whitelist with special detecting.">')
        self.format = input_vcf.format
        self.source = input_vcf.source
        other = f'##Filter_HotSpotdb_Command="python {os.path.abspath(__file__)} --input_vcf {input_vcf.path} --out_vcf {self.path} --hot_db {self.db}'
        self.other = input_vcf.other
        self.other.append(other)
        self.headerline = input_vcf.headerline
        self.contig = input_vcf.contig
        self.mut_df = input_vcf.mut_df
        self.tumor = input_vcf.tumor
        self.normal = input_vcf.normal

    def get_hotspot(self):
        hot_df = pd.read_excel(self.db, names=['基因', '转录本', '位置', '变异', '变异类型', '是否为超级白名单', '备注'])
        hot_df = hot_df[hot_df['是否为超级白名单'] == '是']
        hot_df['转录本'] = hot_df['转录本'].str.split('.', expand=True)[0]
        hot_df['chr'] = hot_df['位置'].str.split(':', expand=True)[0]
        hot_df['start'] = hot_df['位置'].str.split(':', expand=True)[1].str.split('-', expand=True)[0]
        hot_df['end'] = hot_df['位置'].str.split(':', expand=True)[1].str.split('-', expand=True)[1]
        hot_df['AAChange'] = hot_df['基因'].str.cat(hot_df['变异'], sep=':')

        def filterdb(x):
            filter_tags = set(x['filter'].split(';'))
            chrom = x['chr']
            pos = x['pos']
            vaf = float(re.search('AF=(.*?);', x['info']).group(1))
            reads = int(re.search('VD=(.*?);', x['info']).group(1))
            normal_reads = int(x[self.normal].split(':')[2])
            msirep = float(re.search('MSI=(.*?);', x['info']).group(1))
            base_filter = vaf >= 0.004 and reads >= 4
            special_filter = vaf >= 0.004 and reads >= 4 and normal_reads == 0 and not ((msirep >= 3 and vaf <= 0.05) or (msirep >= 5 and vaf <= 0.1))
            gene = re.search('Gene.refGene=(.*?);', x['info']).group(1)
            function = re.search('Func.refGene=(.*?);', x['info']).group(1)
            exon_function = re.search('ExonicFunc.refGene=(.*?);', x['info']).group(1)
            clinvar = re.search('CLNSIG=(.*?);', x['info']).group(1)
            aa_change = re.search('AAChange.refGene=(.*?);', x['info']).group(1)
            # i.split(":")[4][2:]是蛋白质改变（去掉前面的p.与白名单进行匹配）
            mut_set = {f'{i.split(":")[4][2:]}' for i in aa_change.split(',') if 'p.' in i}

            # 遍历超级白名单列表（此功能本质上就是嵌套循环，只能通过itertuples加速）。
            for row in hot_df.itertuples():
                # 进行特殊变异检测匹配（特定区间内的某类突变都算）
                if (  # 先匹配基因与区间
                        getattr(row, '变异类型') == '特殊检测' and getattr(row, '基因') == gene and
                        getattr(row, 'chr') == chrom and getattr(row, 'start') <= pos <= getattr(row, 'end')
                ) and (  # 再分别匹配各变异类型
                        (getattr(row, '变异') == 'fs' and exon_function.startswith('frameshift')) or
                        (getattr(row, '变异') == 'del' and exon_function == 'frameshift_deletion') or
                        (getattr(row, '变异') == 'ins' and exon_function == 'frameshift_insertion') or
                        (getattr(row, '变异') == 'X' and exon_function == 'stopgain') or
                        (getattr(row, '变异') == 'splice' and function == 'splicing') or
                        (getattr(row, '变异') == 'Oncogenic Mutations' and clinvar == 'Pathogenic') or
                        (getattr(row, '变异') == 'Promoter')
                ) and special_filter:  # 特殊变异检测须进行更严格的过滤，尤其是STR导致的fs过滤，否则报出假阳过多
                    filter_tags.add('PASS')
                    if 'germline_risk' in filter_tags:
                        filter_tags.remove('germline_risk')
                    if 'SuperHotSpot' not in x['info']:
                        x['info'] = f"{x['info']};SuperHotSpot_Special"
                # 进行非特殊变异白名单匹配
                elif getattr(row, '基因') == gene and getattr(row, '变异') in mut_set and base_filter:
                    filter_tags.add('PASS')
                    if 'germline_risk' in filter_tags:
                        filter_tags.remove('germline_risk')
                    if 'SuperHotSpot' not in x['info']:
                        x['info'] = f"{x['info']};SuperHotSpot"

            x['filter'] = ';'.join(sorted(list(filter_tags)))
            return x

        self.mut_df = self.mut_df.apply(filterdb, axis=1)
        self.mut_df = self.mut_df[self.mut_df['info'].str.contains('SuperHotSpot') == True]


if __name__ == "__main__":
    args = get_args()
    input_vcf = Vcf('Vardict', args.input_vcf)
    input_vcf.get_header()
    input_vcf.get_mutdf()
    output_vcf = HotSpotVcf('HotSpot', args.output_vcf, input_vcf, args.hot_db)
    output_vcf.get_hotspot()
    output_vcf.write_vcf()
