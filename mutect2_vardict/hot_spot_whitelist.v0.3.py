#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Vardict Vcf文件筛选白名单热点突变，兼容特殊变异检测；支持室间质评降低阈值检测。
v0.3更新：普通(ordinary)白名单十分可靠的也保留。
v0.3.1更新：超级白名单点突变应用碱基质量>=35，以排除包含低碱基质量的假阳
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
    parser.add_argument('--sjzp', '-s', help='ctDNA室间质评超低频突变检测模式', action='store_true')
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

    def get_hotspot(self, sjzp=False):
        hot_df = pd.read_excel(self.db, names=['基因', '转录本', '位置', '变异', '变异类型', '是否超级白名单', '备注', '临床意义'])
        # hot_df = hot_df[hot_df['是否超级白名单'] == '是']
        hot_df['转录本'] = hot_df['转录本'].str.split('.', expand=True)[0]
        hot_df['chr'] = hot_df['位置'].str.split(':', expand=True)[0]
        hot_df['start'] = hot_df['位置'].str.split(':', expand=True)[1].str.split('-', expand=True)[0]
        hot_df['end'] = hot_df['位置'].str.split(':', expand=True)[1].str.split('-', expand=True)[1]
        hot_df['end'].fillna(hot_df['start'], inplace=True)     # start和end相同时，表里可能没有end
        hot_df['AAChange'] = hot_df['基因'].str.cat(hot_df['变异'], sep=':')

        def filterdb(x):
            filter_tags = set(x['filter'].split(';'))
            # 获得变异基本信息
            chrom = x['chr']
            pos = x['pos']
            var_type = re.search('TYPE=(.*?);', x['info']).group(1)
            vaf = float(re.search('AF=(.*?);', x['info']).group(1))
            reads = int(re.search('VD=(.*?);', x['info']).group(1))
            bq = float(x[self.tumor].split(':')[10])                # base quality
            normal_reads = int(x[self.normal].split(':')[2])
            normal_total_reads = int(x[self.normal].split(':')[1])
            normal_vaf = normal_reads / normal_total_reads if normal_total_reads > 0 else 0
            msirep = float(re.search('MSI=(.*?);', x['info']).group(1))
            long_indel = len(x['ref']) > 50 or len(x['alt']) > 50
            vaf_ratio = float(vaf / normal_vaf) > 2 if normal_vaf > 0 else True

            # 获得变异注释信息
            gene = re.search('Gene.refGene=(.*?);', x['info']).group(1)
            function = re.search('Func.refGene=(.*?);', x['info']).group(1)
            exon_function = re.search('ExonicFunc.refGene=(.*?);', x['info']).group(1)
            clinvar = re.search('CLNSIG=(.*?);', x['info']).group(1)
            aa_change = re.search('AAChange.refGene=(.*?);', x['info']).group(1)

            # 设置阈值：超级白名单、普通白名单、是否室间质评
            # 特殊变异检测须进行更严格的过滤，尤其是STR导致的fs过滤，否则报出假阳过多
            if not sjzp:
                # 超级白名单特殊检测阈值
                special_filter = (0.004 <= vaf < 0.9 and reads >= 5 and normal_reads == 0 and bq >= 30 and
                                  not long_indel and
                                  not (msirep >= 3 and vaf <= 0.1))
                # 超级白名单阈值
                if gene == 'BRCA1' or gene == 'BRCA2':
                    base_filter = special_filter    # 如果基因是BRCA1或BRCA2，因超级白名单有大量疑似未经验证位点，故提高阈值
                else:
                    # 超级白名单点突变有类似3个38带一个12的，因此也需bq过滤
                    # 有germline_risk假阳性出现，应加入vaf_ratio
                    base_filter = vaf >= 0.004 and reads >= 4 and bq >= 30 and vaf_ratio
                # 普通白名单特殊检测阈值
                ord_special_filter = (var_type == 'SNV' and vaf >= 0.01 and reads >= 8 and bq >= 30 and normal_reads == 0
                                      and not long_indel
                                      and not (msirep >= 3 and vaf <= 0.1))
                # 普通白名单阈值
                ord_base_filter = var_type == 'SNV' and vaf >= 0.01 and reads >= 8 and bq >= 30 and normal_reads == 0
            else:
                base_filter = vaf >= 0.001 and reads >= 2 and bq >= 30 and vaf_ratio
                ord_base_filter = var_type == 'SNV' and vaf >= 0.005 and reads >= 6 and bq >= 30 and normal_reads == 0
                special_filter = (0.001 <= vaf <= 0.9 and reads >= 2 and normal_reads == 0 and bq >= 30 and
                                  not (msirep >= 3 and vaf <= 0.1))
                ord_special_filter = (var_type == 'SNV' and vaf >= 0.005 and reads >= 6 and bq >= 30 and normal_reads == 0
                                         and not (msirep >= 3 and vaf <= 0.1))

            # i.split(":")[4][2:]是蛋白质改变（去掉前面的p.与白名单进行匹配）
            mut_set = {f'{i.split(":")[4][2:]}' for i in aa_change.split(',') if 'p.' in i}

            # 遍历白名单列表（此功能本质上就是嵌套循环，只能通过itertuples加速）。
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
                ):
                    if getattr(row, '是否超级白名单') == '是' and special_filter:
                        filter_tags.add('PASS')
                        if 'germline_risk' in filter_tags:
                            filter_tags.remove('germline_risk')
                        if 'SuperHotSpot' not in x['info']:
                            x['info'] = f"{x['info']};SuperHotSpot_Special"
                    elif getattr(row, '是否超级白名单') == '否' and ord_special_filter:
                        filter_tags.add('PASS')
                        # if 'germline_risk' in filter_tags:
                        #     filter_tags.remove('germline_risk')
                        if 'HotSpot' not in x['info']:
                            x['info'] = f"{x['info']};HotSpot_Special"
                # 进行非特殊检测白名单匹配
                elif getattr(row, '基因') == gene and getattr(row, '变异') in mut_set:
                    if getattr(row, '是否超级白名单') == '是' and base_filter:
                        filter_tags.add('PASS')
                        if 'germline_risk' in filter_tags:
                            filter_tags.remove('germline_risk')
                        if 'SuperHotSpot' not in x['info']:
                            x['info'] = f"{x['info']};SuperHotSpot"
                    elif getattr(row, '是否超级白名单') == '否' and ord_base_filter:
                        filter_tags.add('PASS')
                        # if 'germline_risk' in filter_tags:
                        #     filter_tags.remove('germline_risk')
                        if 'HotSpot' not in x['info']:
                            x['info'] = f"{x['info']};HotSpot"

            x['filter'] = ';'.join(sorted(list(filter_tags)))
            return x

        self.mut_df = self.mut_df.apply(filterdb, axis=1)
        self.mut_df = self.mut_df[self.mut_df['info'].str.contains('HotSpot') == True]


if __name__ == "__main__":
    args = get_args()
    input_vcf = Vcf('Vardict', args.input_vcf)
    input_vcf.get_header()
    input_vcf.get_mutdf()
    output_vcf = HotSpotVcf('HotSpot', args.output_vcf, input_vcf, args.hot_db)
    output_vcf.get_hotspot(args.sjzp)
    output_vcf.write_vcf()
