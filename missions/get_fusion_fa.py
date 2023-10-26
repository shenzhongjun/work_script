#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
从gtf文件中获得两个融合断点的上下游各500bp位置信息，从fa中提取对应序列
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2022-1-13 16:38:28"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import io
import argparse
import subprocess
import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='从gff文件中获得两个融合断点的上下游各500bp位置信息，从fa中提取对应序列')
    parser.add_argument('--fusion1', help='融合位点A，例如"chr2:42522656"', required=True)
    parser.add_argument('--fusion2', help='融合位点B，例如"chr2:29446394"', required=True)
    parser.add_argument('--gene1', help='融合基因A', required=True)
    parser.add_argument('--gene2', help='融合基因B', required=True)
    parser.add_argument('--gtf', help='gtf文件路径', required=True)
    parser.add_argument('--out', help='输出文件路径', required=True)
    return parser.parse_args()


class FusionGene(object):

    def __init__(self, gene, fusion_loc, gtf_df):
        self.gene = gene
        self.chrom, self.fusion_loc = fusion_loc.split(':')
        self.fusion_loc = int(self.fusion_loc)
        gtf_df['gene'] = gtf_df['attributes'].str.split(';', expand=True)[4].str.strip().str.split(' ', expand=True)[1].str.strip('"')
        gtf_df['transcript'] = gtf_df['attributes'].str.split(';', expand=True)[1].str.strip().str.split(' ', expand=True)[1].str.strip('"')
        self.df = gtf_df[gtf_df['gene'] == self.gene]
        self.df = self.df[(self.df['transcript'] == self.get_long_transcript()) & (self.df['feature'] == 'exon')]

    def get_long_transcript(self):
        trans_df = self.df[self.df['feature'] == 'transcript']
        trans_df['length'] = trans_df.apply(lambda x: x['end'] - x['start'], axis=1)
        trans_df = trans_df.sort_values('length', ascending=False)
        return trans_df.iloc[0]['transcript']

    def get_flank_loc(self, flank=500):
        assert self.chrom == self.df['chrom'].drop_duplicates().iloc[0], 'gtf染色体不匹配，请检查输入参数'
        fusion_interval = []
        exon_intervals = []
        for i, position in enumerate(zip(self.df['start'].tolist(), self.df['end'].tolist())):
            start = position[0]
            end = position[1]
            if start <= self.fusion_loc <= end:
                fusion_interval = [start, end]
                a = self.fusion_loc - flank
                b = self.fusion_loc + flank
                if a < start:
                    exon_intervals.append([start, self.fusion_loc])
                else:
                    exon_intervals.append([a, self.fusion_loc])

                return i, fusion_interval


if __name__ == "__main__":
    args = get_args()
    out = subprocess.getoutput(f'grep -E -w "{args.gene1}|{args.gene2}" {args.gtf}')
    gtf_df = pd.read_table(
        io.StringIO(out),
        names=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes'],
        dtype={'start': int, 'end': int})
    fusion1 = FusionGene(args.gene1, args.fusion1, gtf_df)
    fusion2 = FusionGene(args.gene2, args.fusion2, gtf_df)
    fusion1.df.to_csv('f1.txt', sep='\t', index=False)
    fusion2.df.to_csv('f2.txt', sep='\t', index=False)
    print(list(fusion1.get_flank_loc()))
