#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
VCF的构建、合并与过滤
注：开发版，适配成对样call胚系。当前已更新到使用版；2023年5月24日17:27:34
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-11-7 8:48:11"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import re
import pandas as pd
import gzip


class Vcf(object):
    def __init__(self, tag, path):
        self.tag = tag                              # vcf来源
        self.path = path                            # vcf文件路径
        self.compressed = self.check_file_format()  # vcf是否压缩文件
        self.file_format = ''                       # vcf文件格式
        self.headerline = ''                        # vcf表头
        self.change_headerline = False              # 是否替换表头tumor、normal顺序
        self.filter = []                            # filter注释信息
        self.info = []                              # info注释信息
        self.format = []                            # format注释信息
        self.contig = []                            # contig注释信息
        self.source = []                            # source注释信息
        self.other = []                             # other注释信息
        self.tumor = ''                             # vcf表头中tumor样本名
        self.tumor2 = ''                            # WES+Panel组合中Panel样本名
        self.normal = ''                            # vcf表头中normal样本名
        self.vcf_samples = []                       # vcf表头中所有样本名
        self.mut_df = pd.DataFrame()                # vcf内容

    def get_header(self):
        f = gzip.open(self.path, 'rt') if self.compressed else open(self.path)
        self.file_format = f.readline().strip()

        for line in f:
            if line.startswith('#'):
                line = line.strip()
                if line.startswith('##tumor_sample'):
                    self.tumor = line.strip().split('=')[1]
                    if self.tag != 'Haplotype':     # Haplotype涉及到多个T，单独处理
                        self.other.append(line)
                elif line.startswith('##tumor_sample2'):
                    self.tumor2 = line.strip().split('=')[1]
                elif line.startswith('##normal_sample'):
                    self.normal = line.strip().split('=')[1]
                    self.other.append(line)
                elif line.startswith('##FILTER'):
                    self.filter.append(line)
                elif line.startswith('##source'):
                    self.source.append(line)
                elif line.startswith('##INFO'):
                    self.info.append(line)
                elif line.startswith('##FORMAT'):
                    self.format.append(line)
                elif line.startswith('##contig'):
                    self.contig.append(line)
                elif line.startswith('#CHROM'):
                    if self.tag == 'Mutect2':
                        self.headerline = self.check_header(line)
                    elif self.tag == 'Haplotype':
                        self.headerline = line
                        self.vcf_samples = line.split('\t')[9:]
                    else:
                        self.tumor, self.normal = line.split('\t')[9:]
                        self.headerline = line
                else:
                    self.other.append(line)
        f.close()

    def check_header(self, header):
        """Mutect2的VCF可能Normal在前Tumor在后，需调换位置"""
        cols = header.split('\t')[:9]
        sample1, sample2 = header.split('\t')[9:]
        if sample1 == self.tumor and sample2 == self.normal:
            return header
        elif sample1 == self.normal and sample2 == self.tumor:
            self.change_headerline = True
            print(f'======== {self.tag}-VCF中肿瘤和对照样本位置与预期顺序不符，已修改为肿瘤在前对照在后 ======== ')
            return '\t'.join(cols + [self.tumor, self.normal])
        else:
            raise BaseException(f'======== 文件"{self.path}"的样本与给出的样本名不符，请核对！ ======== ')

    def check_file_format(self):
        if self.path.endswith('.vcf'):
            return False
        elif self.path.endswith('.vcf.gz'):
            return True
        else:
            raise IOError(f'vcf输入文件格式有误:{self.path}')

    def get_mutdf(self):
        if self.tag == 'Haplotype':
            names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + self.vcf_samples
        elif self.change_headerline:
            names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.normal, self.tumor]
        else:
            names = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal]

        if self.compressed:
            self.mut_df = pd.read_table(self.path, compression='gzip', comment='#', names=names, dtype='str')  # 全部列均转为str
        else:
            self.mut_df = pd.read_table(self.path, comment='#', names=names, dtype='str')

        if self.change_headerline:
            self.mut_df = self.mut_df.reindex(
                ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', self.tumor, self.normal], axis=1)

    def write_vcf(self):
        assert not self.compressed, '脚本不支持写入gz文件！'
        with open(self.path, 'w') as w:
            w.write(f'{self.file_format}\n')
            w.write('\n'.join(self.source) + '\n')
            w.write('\n'.join(self.info) + '\n')
            w.write('\n'.join(self.filter) + '\n')
            w.write('\n'.join(self.format) + '\n')
            w.write('\n'.join(self.contig) + '\n' if self.contig else '')   # Vardict原始vcf没有contig
            w.write('\n'.join(self.other) + '\n')
            w.write(f'{self.headerline}\n')
        self.mut_df.drop_duplicates(keep='first', inplace=True)
        self.mut_df.to_csv(self.path, mode='a', header=False, sep='\t', index=False)

    # ======================VCF合并相关=======================
    def handle_mutect_half(self, merged_df):
        """VCF合并：Mutect2独有位点单独处理，亦可用于Mutect2+Vardict与Haplotype合并时Mutect2+Vardict独有位点处理"""
        mutect_half = merged_df.loc[
            merged_df['_merge'] == 'left_only',
            ['chr', 'pos', 'id_x', 'ref', 'alt', 'qual_x', 'filter_x', 'info_x', 'format_x', f'{self.tumor}_x', f'{self.normal}_x']
        ]
        mutect_half = mutect_half.rename(
            columns={'id_x': 'id', 'qual_x': 'qual', 'filter_x': 'filter', 'info_x': 'info', 'format_x': 'format',
                     f'{self.tumor}_x': self.tumor, f'{self.normal}_x': self.normal}
        )
        mutect_half['info'] = mutect_half['info'].apply(lambda x: x + ';SOFT=Mutect2' if 'SOFT=' not in x else x)

        return mutect_half

    def handle_vardict_half(self, merged_df):
        """VCF合并：Vardict独有位点单独处理"""
        vardict_half = merged_df.loc[
            merged_df['_merge'] == 'right_only',
            ['chr', 'pos', 'id_y', 'ref', 'alt', 'qual_y', 'filter_y', 'info_y', 'format_y', f'{self.tumor}_y', f'{self.normal}_y']
        ]
        vardict_half = vardict_half.rename(
            columns={'id_y': 'id', 'qual_y': 'qual', 'filter_y': 'filter', 'info_y': 'info', 'format_y': 'format',
                     f'{self.tumor}_y': self.tumor, f'{self.normal}_y': self.normal}
        )
        vardict_half['info'] += ';SOFT=Vardict'
        return vardict_half

    def handle_haplo_half(self, merged_df):
        """VCF合并：Haplotype独有位点单独处理"""
        haplo_half = merged_df.loc[
            merged_df['_merge'] == 'right_only',
            ['chr', 'pos', 'id_y', 'ref', 'alt', 'qual_y', 'filter_y', 'info_y', 'format_y', f'{self.tumor}_y', f'{self.normal}_y']
        ]
        haplo_half = haplo_half.rename(
            columns={'id_y': 'id', 'qual_y': 'qual', 'filter_y': 'filter', 'info_y': 'info', 'format_y': 'format',
                     f'{self.tumor}_y': self.tumor, f'{self.normal}_y': self.normal}
        )
        haplo_half['info'] += ';SOFT=Haplotype'
        haplo_half['id'] = '.'  # rsIDs可以通过注释软件获得，且有rsID干扰vep的'Uploaded_variation'列
        haplo_half[self.tumor] = haplo_half[self.tumor].apply(lambda x: '1/1:0,0:0:0:0' if './.' in x else x)
        return haplo_half

    @staticmethod
    def merge_mh(x, t, n):
        """VCF合并：Mutect2(+Vardict)和Haplotype交集位点合并，以Mutect2格式为主"""
        x['id'] = '.'  # rsIDs可以通过注释软件获得
        x['qual'] = x['qual_y']
        x['filter'] = x['filter_x']
        soft = 'SOFT=Haplotype,Mutect2' if 'SuperHotSpot' not in x['info_x'] else 'SOFT=Haplotype,Mutect2,Vardict'
        x['info'] = f"{x['info_x']};{soft}" if 'SuperHotSpot' not in x['info_x'] else re.sub('SOFT=.*', soft, x['info_x'])
        x['format'] = x['format_x']
        # x[t] = f'.:{x[f"{t}_x"].split(":")[1]}:{x[f"{t}_x"].split(":")[3]}:.:.'   # 以Haplotype为主时用Mutect2结果填充t_format
        x[t] = x[f'{t}_x']
        x[n] = x[f'{n}_x']
        return x

    # Mutect2和Vardict交集位点合并
    @staticmethod
    def merge_mv(x, t, n):
        """VCF合并：Mutect2和Vardict交集位点合并，以Mutect2为主"""
        x['id'] = '.'
        x['qual'] = x['qual_y']  # Mutect2如果不做VQSR则没有Q值
        x['filter'] = ';'.join(
            sorted(list(set(x['filter_x'].split(';')).union(set(x['filter_y'].split(';'))))))  # filter直接合并
        x['info'] = x['info_x'] + ';' + re.sub('DP=.*?;', '', x['info_y']) + ';SOFT=Mutect2,Vardict'  # info去掉Vardict的DP后其余信息放在Mutect2后面
        x['format'] = x['format_x']
        x[t] = x[f'{t}_x']
        x[n] = x[f'{n}_x']
        return x

    # --------VCF过滤相关--------
    @staticmethod
    def check_var_origin(x, n):
        """合并后VCF确定突变来源"""
        filter_tags = set(x['filter'].split(';'))
        soft = re.search('SOFT=(.*)', x['info']).group(1)
        if x['format'].startswith('GT:AD:DP'):      # Haplotype格式FORMAT列
            vaf = float(int(x[n].split(':')[1].split(',')[1])/int(x[n].split(':')[2]))
        else:       # Mutect2格式FORMAT列
            vaf = float(x[n].split(':')[2])

        if 'Haplotype' in soft and vaf >= 0.25:
            filter_tags.add('germline')
            filter_tags.add('PASS')
            if 'germline_risk' in filter_tags:
                filter_tags.remove('germline_risk')
        elif soft == 'Haplotype' and vaf < 0.25:
            filter_tags.add('germline_risk')
        elif soft == 'Haplotype,Mutect2' and vaf <= 0.25 and 'germline_risk' in filter_tags:
            filter_tags.add('germline_risk')

        x['filter'] = ';'.join(sorted(list(filter_tags)))
        return x

