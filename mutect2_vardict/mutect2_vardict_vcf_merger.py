#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mutect2和Vardict过滤后vcf结果合并，输出vcf文件。
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-11-7 8:48:11"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import re
import argparse
from handle_vcf import Vcf


def get_args():
    parser = argparse.ArgumentParser(description='Mutect2和Vardict过滤后vcf结果合并，输出vcf文件。')
    parser.add_argument('--mutect2_vcf', '-m', help='Mutect2 *.vcf.gz输入文件路径', required=True)
    parser.add_argument('--vardict_vcf', '-v', help='Vardict *.vcf.gz输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf合并结果输出文件路径', required=True)
    return parser.parse_args()


class Merger(Vcf):
    def __init__(self, tag, path, mutect_vcf, vardict_vcf):
        super().__init__(tag, path)
        self.mutect_vcf = mutect_vcf
        self.vardict_vcf = vardict_vcf
        assert self.mutect_vcf.tumor == self.vardict_vcf.tumor, 'Mutect2-VCF和Vardict-VCF的肿瘤样本不一致！'
        assert self.mutect_vcf.tumor == self.vardict_vcf.tumor, 'Mutect2-VCF和Vardict-VCF的对照样本不一致！'
        self.tumor = self.vardict_vcf.tumor
        self.normal = self.vardict_vcf.normal

    def merge(self):
        self.file_format = self.merge_file_format()
        self.filter = self.merge_vcf_filter()
        self.info = self.merge_vcf_info()
        self.format = self.merge_vcf_format()
        self.source = self.merge_vcf_source()
        # self.contig = self.merge_vcf_contig()
        self.other = self.merge_vcf_other()
        self.headerline = self.merge_headerline()
        self.contig = self.mutect_vcf.contig
        self.mut_df = self.merge_mutdf()
        self.write_vcf()

    def merge_file_format(self):
        """合并vcf表头之file_format"""
        if self.mutect_vcf.file_format == self.vardict_vcf.file_format:
            return self.mutect_vcf.file_format
        else:
            raise BaseException('======== 两个vcf文件版本不一致，请确认是否来自最新版Mutect2和Vardict！======== ')

    def merge_headerline(self):
        """合并vcf表头之headerline"""
        if self.mutect_vcf.headerline == self.vardict_vcf.headerline:
            return self.mutect_vcf.headerline
        else:
            raise BaseException('======== 两个vcf文件格式或样本不一致，请核对！======== ')

    def merge_vcf_filter(self):
        """合并vcf表头之filter"""
        merged_filer = list(set(self.mutect_vcf.filter).union(set(self.vardict_vcf.filter)))
        merged_filer.sort()
        return merged_filer

    def merge_vcf_info(self):
        """合并vcf表头之info，仅有DP重复，以Vardict为主"""
        mutect_info = [i for i in self.mutect_vcf.info if i.split(',')[0].split('=')[-1] != 'DP']
        mutect_info.append('##INFO=<ID=SOFT,Number=1,Type=String,Description="Which software called this variant">')
        return self.vardict_vcf.info + mutect_info

    def merge_vcf_format(self):
        """合并vcf表头之format，AD、GT、DP、AF重复，以Vardict为主"""
        mutect_format = [i for i in self.mutect_vcf.format if i.split(',')[0].split('=')[-1] not in ['AD', 'GT', 'DP', 'AF']]
        return self.vardict_vcf.format + mutect_format

    def merge_vcf_source(self):
        """合并vcf表头之source"""
        return self.mutect_vcf.source + self.vardict_vcf.source

    def merge_vcf_contig(self):
        """合并vcf表头之contig，以Mutect2为主不合并也可以"""
        return self.mutect_vcf.contig + self.vardict_vcf.contig

    def merge_vcf_other(self):
        """合并vcf表头之other"""
        info = f'##Mutect2_Vardict_vcf_mergeCommand="python {os.path.abspath(__file__)} --mutect2_vcf {self.mutect_vcf.path}' \
               f' --vardict_vcf {self.vardict_vcf.path} --output_vcf {self.path}'
        self.vardict_vcf.other.append(info)
        return self.mutect_vcf.other + self.vardict_vcf.other

    def merge_mutdf(self):
        """合并vcf突变"""
        merged_df = self.mutect_vcf.mut_df.merge(
            self.vardict_vcf.mut_df,
            how='outer',
            on=['chr', 'pos', 'ref', 'alt'],
            sort=True,
            indicator=True
        )
        # 为节省资源，取Mutect2独有位点单独处理
        merged_df_left = merged_df.loc[
            merged_df['_merge'] == 'left_only',
            ['chr', 'pos', 'id_x', 'ref', 'alt', 'qual_x', 'filter_x', 'info_x', 'format_x',
             f'{self.tumor}_x', f'{self.normal}_x']
        ]
        merged_df_left = merged_df_left.rename(
            columns={'id_x': 'id', 'qual_x': 'qual', 'filter_x': 'filter', 'info_x': 'info', 'format_x': 'format',
                     f'{self.tumor}_x': self.tumor, f'{self.normal}_x': self.normal}
        )
        merged_df_left['info'] += ';SOFT=Mutect2'

        # 为节省资源，取Vardict独有位点单独处理
        merged_df_right = merged_df.loc[
            merged_df['_merge'] == 'right_only',
            ['chr', 'pos', 'id_y', 'ref', 'alt', 'qual_y', 'filter_y', 'info_y', 'format_y',
             f'{self.tumor}_y', f'{self.normal}_y']
        ]
        merged_df_right = merged_df_right.rename(
            columns={'id_y': 'id', 'qual_y': 'qual', 'filter_y': 'filter', 'info_y': 'info', 'format_y': 'format',
                     f'{self.tumor}_y': self.tumor, f'{self.normal}_y': self.normal}
        )
        merged_df_right['info'] += ';SOFT=Vardict'

        # 取两款软件交集位点进行合并
        merged_df_both = merged_df.loc[merged_df['_merge'] == 'both']
        merged_df_both = merged_df_both.apply(self.mergex, axis=1, args=(self.tumor, self.normal))
        merged_df_both = merged_df_both[['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]]
        merged_df_both['info'] += ';SOFT=Vardict,Mutect2'

        # 独有位点与交集位点合并后列顺序重新排列
        merged_df_result = merged_df_right.append([merged_df_both, merged_df_left], ignore_index=True, sort=True)
        merged_df_result = merged_df_result.reindex(
            ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal],
            axis=1
        )
        return merged_df_result

    @staticmethod
    def mergex(x, t, n):
        x['id'] = '.'   # rsIDs可以通过注释软件获得，Mutect2和Vardict结果中均无rsID
        x['qual'] = x['qual_y']     # Mutect2如果不做VQSR则没有Q值
        # x['filter'] = ';'.join(sorted(list(set(x['filter_x'].split(';')).union(set(x['filter_y'].split(';'))))))    # filter直接合并
        x['filter'] = x['filter_x'] if 'germline' not in x['filter_y'] else x['filter_y']
        x['info'] = x['info_y'] + ';' + re.sub('DP=.*?;', '', x['info_x'])      # info去掉Mutect2的DP后其余信息放在Vardict后面
        x['format'] = x['format_y'] + re.sub('GT:AD:AF:DP', '', x['format_x'])         # format去掉Mutect2的AD、GT、DP、AF后其余信息放在Vardict后面
        x[t] = x[f'{t}_y'] + ':' + ':'.join(x[f'{t}_x'].split(':')[4:])
        x[n] = x[f'{n}_y'] + ':' + ':'.join(x[f'{n}_x'].split(':')[4:])
        return x


if __name__ == "__main__":
    args = get_args()

    mutect2_vcf = Vcf('Mutect2', args.mutect2_vcf)
    mutect2_vcf.get_header()
    mutect2_vcf.get_mutdf()

    vardict_vcf = Vcf('Vardict', args.vardict_vcf)
    vardict_vcf.get_header()
    vardict_vcf.get_mutdf()

    mergedvcf = Merger('Merger', args.output_vcf, mutect2_vcf, vardict_vcf)
    mergedvcf.merge()
