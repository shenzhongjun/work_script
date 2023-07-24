#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Mutect2 vcf、Vardict超级白名单vcf和Haplotype成对样胚系变异vcf结果合并，输出vcf文件。
v0.3更新：增加对WES+Panel双tumor样本call胚系变异后合并的支持
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-2-8"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import re
import argparse
import pandas as pd
from handle_vcf_2 import Vcf


def get_args():
    parser = argparse.ArgumentParser(description='Mutect2 vcf、Vardict超级白名单vcf和Haplotype胚系变异vcf结果合并，输出vcf文件。')
    parser.add_argument('--mutect2_vcf', '-m', help='Mutect2 vcf输入文件路径', required=True)
    parser.add_argument('--vardict_vcf', '-v', help='Vardict vcf输入文件路径', required=True)
    parser.add_argument('--haplo_vcf', '-t', help='Haplotype vcf输入文件路径', required=True)
    parser.add_argument('--output_vcf', '-o', help='vcf合并结果输出文件路径', required=True)
    return parser.parse_args()


class Merger(Vcf):
    def __init__(self, tag, path, mutect_vcf, vardict_vcf, haplo_vcf):
        super().__init__(tag, path)
        self.mutect_vcf = mutect_vcf
        self.vardict_vcf = vardict_vcf
        self.haplo_vcf = haplo_vcf
        self.nohotspot = True if self.vardict_vcf.mut_df.empty else False
        self.tumor = self.mutect_vcf.tumor
        self.normal = self.mutect_vcf.normal
        if self.haplo_vcf.tumor2:   # 如果使用WES+Panel两个T样本参与胚系变异检出
            if self.tumor == self.haplo_vcf.tumor:  # Mutect2的T和Haplotype的T相同，是WES样本
                self.haplo_vcf.mut_df = self.haplo_vcf.mut_df.drop(self.haplo_vcf.tumor2, axis=1)
            elif self.tumor == self.haplo_vcf.tumor2:   # Mutect2的T和Haplotype的T2相同，是Panel样本
                self.haplo_vcf.mut_df = self.haplo_vcf.mut_df.drop(self.haplo_vcf.tumor, axis=1)
            else:
                raise BaseException('！！！Mutect2肿瘤样本和Haplotype肿瘤样本不对应！！！')
            self.haplo_vcf.mut_df = self.haplo_vcf.mut_df.reindex(
                ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal], axis=1
            )

        assert self.mutect_vcf.tumor == self.vardict_vcf.tumor, 'Mutect2和Vardict的肿瘤样本不一致！'
        assert self.mutect_vcf.normal == self.vardict_vcf.normal == self.haplo_vcf.normal, 'Mutect2、Vardict、Haplotype的对照样本不一致！'


    def merge(self):
        self.file_format = self.merge_file_format()
        self.filter = self.merge_vcf_filter()
        self.info = self.merge_vcf_info()
        self.format = self.merge_vcf_format()
        self.source = self.merge_vcf_source()
        self.contig = self.mutect_vcf.contig
        self.other = self.merge_vcf_other()
        self.headerline = self.mutect_vcf.headerline
        self.mut_df = self.merge_mutdf()

    def merge_file_format(self):
        """合并vcf表头之file_format"""
        if self.mutect_vcf.file_format == self.vardict_vcf.file_format == self.haplo_vcf.file_format:
            return self.mutect_vcf.file_format
        else:
            raise BaseException('！！！Vcf文件版本不一致，请使用最新版GATK和Vardict！！！')

    def merge_vcf_filter(self):
        """合并vcf表头之filter"""
        if self.nohotspot:
            merged_filer = list(set(self.mutect_vcf.filter).union(set(self.haplo_vcf.filter)))
        else:
            merged_filer = list(set(self.mutect_vcf.filter).union(set(self.vardict_vcf.filter)).union(set(self.haplo_vcf.filter)))
        merged_filer.append('##FILTER=<ID=germline,Description="Evidence show this site is germline">')
        merged_filer.sort()
        return merged_filer

    def merge_vcf_info(self):
        """合并vcf表头之info，DP和SOR有重复，以GATK为主"""
        if self.nohotspot:
            merged_info = self.mutect_vcf.info + self.haplo_vcf.info
        else:
            vardict_info = [i for i in self.vardict_vcf.info if i.split(',')[0].split('=')[-1] not in ['DP', 'SOR']]
            merged_info = vardict_info + self.mutect_vcf.info + self.haplo_vcf.info
        merged_info.append('##INFO=<ID=SOFT,Number=1,Type=String,Description="Which software called this variant">')
        return merged_info

    def merge_vcf_format(self):
        """合并vcf表头之format，AD、GT、DP、AF重复，以Mutect2为主，Haplotype与Mutect2相同"""
        if self.nohotspot:
            return self.mutect_vcf.format
        else:
            vardict_format = [i for i in self.vardict_vcf.format if i.split(',')[0].split('=')[-1] not in ['AD', 'GT', 'DP', 'AF']]
            return self.mutect_vcf.format + vardict_format

    def merge_vcf_source(self):
        """合并vcf表头之source"""
        if self.nohotspot:
            return self.mutect_vcf.source + self.haplo_vcf.source
        else:
            return self.mutect_vcf.source + self.vardict_vcf.source + self.haplo_vcf.source

    def merge_vcf_other(self):
        """合并vcf表头之other"""
        info = f'##Mutect2_Vardict_Haplotype_vcf_mergeCommand="python {os.path.abspath(__file__)} --mutect2_vcf {self.mutect_vcf.path}' \
               f' --vardict_vcf {self.vardict_vcf.path} --haplo_vcf {self.haplo_vcf.path} --output_vcf {self.path}'
        self.haplo_vcf.other.append(info)
        if self.nohotspot:
            return self.mutect_vcf.other + self.haplo_vcf.other
        else:
            return self.mutect_vcf.other + self.vardict_vcf.other + self.haplo_vcf.other

    def merge_mutdf(self):
        """合并vcf突变并判断突变来源"""
        if self.nohotspot:
            # 无Vardict结果，Mutect2体细胞突变位点与Haplotype胚系突变位点进行合并
            merged_df = self.mutect_vcf.mut_df.merge(self.haplo_vcf.mut_df, how='outer', on=['chr', 'pos', 'ref', 'alt'],
                                                     sort=True, indicator=True)

            # 取Mutect2独有位点和Haplotype独有位点
            mutect_half = self.handle_mutect_half(merged_df)
            haplo_half = self.handle_haplo_half(merged_df)

            # Mutect2和Haplotype交集位点，以Mutect2格式为主进行合并
            merged_df_both = merged_df.loc[merged_df['_merge'] == 'both']
            if not merged_df_both.empty:
                merged_df_both = merged_df_both.apply(Vcf.merge_mh, axis=1, args=(self.tumor, self.normal))
                merged_df_both = merged_df_both[
                    ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]]
            else:
                merged_df_both = pd.DataFrame(
                    columns=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]
                )

            # 独有位点与交集位点合并后列顺序重新排列
            merged_df_result = mutect_half.append([merged_df_both, haplo_half], ignore_index=True, sort=True)
            merged_df_result = merged_df_result.reindex(
                ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal], axis=1
            )
        else:
            # 有Vardict结果，Mutect2、Vardict体细胞突变位点与Haplotype胚系突变位点进行合并
            tmp_df = self.mutect_vcf.mut_df.merge(self.vardict_vcf.mut_df, how='outer', on=['chr', 'pos', 'ref', 'alt'],
                                                  sort=True, indicator=True)
            # 取Mutect2独有位点和Vardict独有位点
            mutect_half = self.handle_mutect_half(tmp_df)
            vardict_half = self.handle_vardict_half(tmp_df)

            # Mutect2和Vardict交集位点进行合并
            tmp_df_both = tmp_df.loc[tmp_df['_merge'] == 'both']
            if not tmp_df_both.empty:
                tmp_df_both = tmp_df_both.apply(Vcf.merge_mv, axis=1, args=(self.tumor, self.normal))
                tmp_df_both = tmp_df_both[
                    ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]]
            else:
                tmp_df_both = pd.DataFrame(
                    columns=['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]
                )

            # Mutect2和Vardict独有位点与交集位点合并后列顺序重新排列
            tmp_df_result = mutect_half.append([tmp_df_both, vardict_half], ignore_index=True, sort=True)
            tmp_df_result = tmp_df_result.reindex(
                ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal],
                axis=1
            )

            # Mutect2和Vardict合并结果与Haplotype结果合并
            merged_df = tmp_df_result.merge(self.haplo_vcf.mut_df, how='outer', on=['chr', 'pos', 'ref', 'alt'],
                                            sort=True, indicator=True)

            # 取Mutect2+Vardict独有位点和Haplotype独有位点
            mut_var_half = self.handle_mutect_half(merged_df)
            haplo_half = self.handle_haplo_half(merged_df)

            # Mutect2+Vardict和Haplotype交集位点，以Haplotype格式为主进行合并
            merged_df_both = merged_df.loc[merged_df['_merge'] == 'both']
            merged_df_both = merged_df_both.apply(Vcf.merge_mh, axis=1, args=(self.tumor, self.normal))
            if not merged_df_both.empty:        # 小芯片可能没有交集
                merged_df_both = merged_df_both[
                    ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]]

            # 独有位点与交集位点合并后列顺序重新排列
            merged_df_result = mut_var_half.append([merged_df_both, haplo_half], ignore_index=True, sort=True)
            merged_df_result = merged_df_result.reindex(
                ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal], axis=1
            )

        merged_df_result = merged_df_result.apply(Vcf.check_var_origin, axis=1, args=(self.normal, ))
        return merged_df_result


if __name__ == "__main__":
    args = get_args()

    mutect2_vcf = Vcf('Mutect2', args.mutect2_vcf)
    mutect2_vcf.get_header()
    mutect2_vcf.get_mutdf()

    vardict_vcf = Vcf('Vardict', args.vardict_vcf)
    vardict_vcf.get_header()
    vardict_vcf.get_mutdf()

    haplo_vcf = Vcf('Haplotype', args.haplo_vcf)
    haplo_vcf.get_header()
    haplo_vcf.get_mutdf()

    mergedvcf = Merger('Merger', args.output_vcf, mutect2_vcf, vardict_vcf, haplo_vcf)
    mergedvcf.merge()
    mergedvcf.write_vcf()
