#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量从vcf格式变异结果中提取阳性位点Call变异情况，需要给出SAMPLENAME.xxx.vcf(.gz)格式的文件
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2021-10-29 10:16:38"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import glob
import gzip
import os
import re
from io import StringIO

import numpy as np
import pandas as pd

pd.reset_option('display.float_format')


def get_args():
    parser = argparse.ArgumentParser(description='批量从vcf格式或annovar注释结果中提取阳性位点Call变异情况')
    parser.add_argument('--soft', '-f',
                        help='vcf文件来源，默认merged（vardict+mutect2），从sentieon、mutect2、vardict、merged中选择',
                        choices=['merged', 'vardict', 'sentieon', 'mutect2'], default='merged')
    parser.add_argument('--glob_path', '-g', help='vcf文件路径，支持"*"')
    parser.add_argument('--sample_list', '-s', help='vcf文件路径，以sample_list.txt提供，第一列为sample、第二列为path')
    parser.add_argument('--positive_list', '-p',
                        help='阳性位点列表文件，类似vcf格式：chr	start	end	ref	alt	gene',
                        required=True)
    parser.add_argument('--out', '-o', help='输出文件路径', required=True)
    parser.add_argument('--annoed', '-a', help='输入文件是UMI-vcf注释结果', action='store_true')
    parser.add_argument('--mrd', '-m', help='输出MRD结果，af保留六位小数', action='store_true')
    parser.add_argument('--mrd_product', '-mp', help='输出MRD产品结果，非PASS位点af设为0，无样本行', action='store_true')
    return parser.parse_args()


def get_tumor_name(pathname):
    file_prefix = os.path.split(pathname)[1].split('.')[0]
    if file_prefix == 'snv_indel':  # 研发流程出的结果文件不包含样本名
        f = gzip.open(pathname, 'rt') if 'gz' in pathname else open(pathname, 'rt')
        for line in f:
            if line.startswith('##tumor_sample'):
                return line.strip().split('=')[1]
        f.close()
    else:
        return file_prefix


def get_file_paths(pathname):
    paths = glob.glob(pathname)
    p_dic = {}
    for full_path in paths:
        library_n = get_tumor_name(full_path)
        p_dic[library_n] = full_path
    return p_dic


def get_file_paths2(sample_list):
    with open(sample_list) as f:
        p_dic = {}
        for line in f:
            library_n, full_path = line.strip().split('\t')
            p_dic[library_n] = full_path
        return p_dic


def get_umi_library_df(file_path, library_n):
    df = pd.read_table(file_path, dtype={'Start': 'str'})
    df['Position'] = df['Chr'].str.cat(df['Start'], sep='_')
    df = df[['Position', 'TotalDepth', 'AltDepth', 'MutFreq']]
    df = df.rename(columns={'TotalDepth': f'TotalDepth_{library_n}'})
    df = df.rename(columns={'TotalDepth': f'总深度_{library_n}', 'AltDepth': '突变深度', 'MutFreq': '突变率'})
    return df


def get_vcf_library_df(file_path, library_n, positive_sites, soft='vardict', check_tri=False, mrd=False, mrd_producrt=False):
    """
    Sentieon VCF格式
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
    ##FORMAT=<ID=AFDP,Number=1,Type=Integer,Description="Read depth to calculate AF">
    vcf FORMAT：
    GT:AD:AF:AFDP:ALTHC:ALT_F1R2:ALT_F2R1:ClippingRankSumPS:DPHC:FOXOG:QSS:REF_F1R2:REF_F2R1:ReadPosEndDistPS:ReadPosRankSumPS
    vcf TUMOR：
    0/1:26151,83:0.003:26234:83:83:0:-3.522:0:26164:1:0.16:36.769:972124,3019:26150:1:36.785:0.818
    check_tri: 如果为True，则合并Position列和Alt列，以排除不需要的triallelic_site

    Vardict VCF格式
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
    ##FORMAT=<ID=VD,Number=1,Type=Integer,Description="Variant Depth">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=ALD,Number=2,Type=Integer,Description="Variant forward, reverse reads">
    ##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference forward, reverse reads">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
    GT:DP:VD:ALD:RD:AD:AF:BIAS:PMEAN:PSTD:QUAL:QSTD:SBF:ODDRATIO:MQ:SN:HIAF:ADJAF:NM
    1/0:120:63:28,35:33,24:57,63:0.525:2,2:50.2:1:37.4:1:0.14931:1.71095:48:126:0.525:0:1.4

    Mutect2 VCF格式
    GT:AD:AF:DP:F1R2:F2R1:FAD:SB
    0/1:275,14:0.093:289:36,23:82,29:173,86:128,147,7,7
    """
    tmp2 = StringIO()
    if 'gz' in file_path:
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp2.write(line)
    else:
        with open(file_path, 'rt') as f:
            for line in f:
                if not line.startswith('##'):
                    tmp2.write(line)
    df = pd.read_table(StringIO(tmp2.getvalue()), dtype={'POS': 'str'})
    result_df = pd.DataFrame(columns=['Position', 'TotalDepth', 'AltDepth', 'MutFreq', 'FILTER'])

    if check_tri:
        df['Position'] = result_df['Position'] = df['#CHROM'].str.cat([df['POS'], df['ALT']], sep='_')
    else:
        df['Position'] = result_df['Position'] = df['#CHROM'].str.cat(df['POS'], sep='_')

    df = df[df['Position'].isin(positive_sites)]

    def get_alt_depth(x):
        if 'F1R2' in x['FORMAT']:
            return int(x[library_n].split(':')[1].split(',')[1])
        else:
            return int(x[library_n].split(':')[5].split(',')[1])

    def get_total_depth(x):
        if 'F1R2' in x['FORMAT']:
            return int(x[library_n].split(':')[3])
        else:
            return int(x[library_n].split(':')[1])

    def get_mutfreq(x):
        return f"{x['AltDepth'] / x['TotalDepth'] * 100:.3f}"

    def get_mrd_product_mutfreq(x):
        if x['FILTER'] == 'PASS':
            return f"{x['AltDepth'] / x['TotalDepth'] * 100:.4f}"
        else:
            return 0

    def get_mrd_mutfreq(x):
        return f"{x['AltDepth'] / x['TotalDepth'] * 100:.4f}"

    def get_filter(x):
        if 'germline_risk' in x['FILTER'] or 'P0.05' in x['FILTER']:
            if 'F1R2' in x['FORMAT']:
                normal_total_depth = int(x.iloc[-2].split(':')[3])
                normal_alt_depth = int(x.iloc[-2].split(':')[1].split(',')[0])
            else:
                normal_total_depth = int(x.iloc[-2].split(':')[1])
                normal_alt_depth = int(x.iloc[-2].split(':')[5].split(',')[1])
            normal_mutfreq = f'{normal_alt_depth / normal_total_depth * 100:.3f}' if normal_total_depth > 0 else 0
            x['FILTER'] = x['FILTER'].replace(
                'germline_risk',
                f'germline_risk({normal_total_depth}/{normal_alt_depth}/{normal_mutfreq})'
            ).replace(
                'P0.05',
                f'P0.05({normal_total_depth}/{normal_alt_depth}/{normal_mutfreq})'
            )
        return x['FILTER']

    if soft == 'merged':
        result_df['AltDepth'] = df.apply(get_alt_depth, axis=1)
        result_df['TotalDepth'] = df.apply(get_total_depth, axis=1)
        result_df['MutFreq'] = result_df.apply(get_mutfreq, axis=1)
    elif soft == 'sentieon':
        result_df['AltDepth'] = df[library_n].str.split(':', expand=True)[1].str.split(',', expand=True)[1].astype(np.int64)
        result_df['TotalDepth'] = df[library_n].str.split(':', expand=True)[3].astype(np.int64)
        result_df['MutFreq'] = result_df.apply(get_mutfreq, axis=1)
    else:
        pass

    result_df['FILTER'] = df.apply(get_filter, axis=1)

    if mrd_producrt:
        result_df['MutFreq'] = result_df.apply(get_mrd_product_mutfreq, axis=1)
    elif mrd:
        result_df['MutFreq'] = result_df.apply(get_mrd_mutfreq, axis=1)

    # result_df['TotalDepth'] = result_df['TotalDepth'].astype('str')
    result_df[['TotalDepth', 'AltDepth', 'MutFreq']] = result_df[['TotalDepth', 'AltDepth', 'MutFreq']].astype('str')
    result_df = result_df.rename(columns={'TotalDepth': f'总深度_{library_n}', 'AltDepth': '突变深度',
                                          'MutFreq': '突变率(%)', 'FILTER': '是否通过过滤'})
    return result_df


def merge_df(positive_df, annoed, mrd, mrd_producrt, soft, glob_path=None, sample_list=None):
    check_tri = False  # 如果该位点处有多等位基因突变，需要指定要查看的突变碱基
    if 'alt' in positive_df.columns:
        positive_df['Position'] = positive_df['Position'].str.cat(positive_df['alt'], sep='_')
        check_tri = True
    if glob_path:
        mut_paths = get_file_paths(glob_path)
    else:
        mut_paths = get_file_paths2(sample_list)

    positive_sites = positive_df['Position'].to_list()
    for library in mut_paths:
        if annoed:
            library_df = get_umi_library_df(mut_paths[library], library)
        else:
            library_df = get_vcf_library_df(mut_paths[library], library, positive_sites, soft, check_tri, mrd, mrd_producrt)
        positive_df = positive_df.merge(library_df, how='left', on='Position')
    return positive_df.drop(columns=['Position']).fillna('-').astype('str')


def write_txt(merged_df, out, mrd_product):
    header = ''
    w = open(out, 'w')
    if not mrd_product:
        for i in merged_df.columns:
            header += '\t' if '总深度' not in i else f'{re.match("总深度_(.*)", i).group(1)}\t'
        w.write(header + '\n')
    header2 = re.sub('总深度_.+?\t', '总深度\t', '\t'.join(merged_df.columns.tolist())). \
        replace('突变深度_x', '突变深度').replace('突变深度_y', '突变深度'). \
        replace('突变率(%)_x', '突变率(%)').replace('突变率(%)_y', '突变率(%)'). \
        replace('是否通过过滤_x', '是否通过过滤').replace('是否通过过滤_y', '是否通过过滤')
    w.write(header2 + '\n')
    for line in merged_df.values.tolist():
        w.write('\t'.join(line) + '\n')
    w.close()


if __name__ == "__main__":
    pd.set_option('max_columns', 15)
    argv = get_args()
    # positive_df header参考：chr	start	end	ref	alt	gene
    positive_df = pd.read_table(argv.positive_list)
    positive_df['Position'] = positive_df[positive_df.columns.values[0]].str.cat(
        positive_df[positive_df.columns.values[1]].astype('str'), sep='_')
    if argv.glob_path:
        merged_df = merge_df(positive_df, argv.annoed, argv.mrd, argv.mrd_product, argv.soft, glob_path=argv.glob_path)
    elif argv.sample_list:
        merged_df = merge_df(positive_df, argv.annoed, argv.mrd, argv.mrd_product, argv.soft,
                             sample_list=argv.sample_list)
    else:
        raise Exception("请提供变异文件路径！")
    write_txt(merged_df, argv.out, argv.mrd_product)
