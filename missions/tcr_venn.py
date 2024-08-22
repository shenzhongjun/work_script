#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
关旭7例样本不同样本类型间TCR克隆Venn图
"""
__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-4-9 09:36:18"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import pandas as pd
from matplotlib import pyplot as plt
import venn

pd.options.display.max_columns = 16

def get_merged_df(sample, sample_type, samples, info):
    df = pd.DataFrame()
    for sample3 in samples:
        for l in info:
            sample4, path, library, _ = l.strip().split('\t')
            if sample4 == sample3:
                df2 = pd.read_table(fr'{path}/{sample4}/05.count_vdj/outs/{sample4}_{library}.merge_vdj_clonetype_scRNACluster.xls')
                df2 = df2.loc[
                    (~df2['aaSeqCDR3_TRA'].isna()) & (~df2['aaSeqCDR3_TRB'].isna()) & (df2['celltype'] == 'T cells')]
                df = df._append(df2)
    if len(df) > 0:
        df['aaseq'] = df['aaSeqCDR3_TRA'] + df['aaSeqCDR3_TRB']
        # duplicates = df.duplicated(subset=['aaseq'], keep=False)  # 标记所有重复值
        # num_duplicates = duplicates.sum()
        # if num_duplicates > 0:
        #     print("存在重复值，共有 {} 个重复值。".format(num_duplicates))
        #     duplicate_rows = df[duplicates]
        #     print("重复值的数据行：")
        #     print(duplicate_rows)
        # else:
        #     print("不存在重复值。")
        df_unique = df.drop_duplicates(subset=['aaseq'], keep='first')
        df_unique.to_csv(fr'venn_plot/{sample}_{sample_type}_merge_vdj_clonetype_scRNACluster.xls', sep='\t', header=True, index=False)
        return set(df['aaseq'].tolist())
    else:
        return set()


sample_list = ["第一例", "第二例", "第三例", "第四例", "第六例", "第七例"]
sample_num = [1, 2, 3, 4, 6, 7]
for sample, number in zip(sample_list, sample_num):
    tumor = []
    normal = []
    ln = []
    bd = []
    with open('samplelist.txt') as f, open('sampleinfo.txt') as info:
        f.readline()
        f = f.readlines()
        info = info.readlines()
        for i in f:  #
            sample2, sample_type, sample_num = i.strip().split('\t')
            # print(sample2, sample_type, sample_num)
            if sample2 == sample:
                if sample_type == '肿瘤T':
                    tumor.append(sample_num)
                elif sample_type == '癌旁N':
                    normal.append(sample_num)
                elif sample_type == '淋巴结LN':
                    ln.append(sample_num)
                elif sample_type == '外周血BD':
                    bd.append(sample_num)
            # print(sample2, tumor, normal, ln, bd)
        tumor_barcodes = get_merged_df(sample, '肿瘤T', tumor, info)
        normal_barcodes = get_merged_df(sample, '癌旁N', normal, info)
        ln_barcodes = get_merged_df(sample, '淋巴结LN', ln, info)
        bd_barcodes = get_merged_df(sample, '外周血BD', bd, info)

        plt.figure(figsize=(6, 6))
        labels = venn.get_labels([tumor_barcodes, normal_barcodes, ln_barcodes, bd_barcodes], fill=['number'])
        fig, ax = venn.venn4(labels, names=['T', 'N', 'LN', 'BD'])
        plt.title(f"Venn Diagram of Sample{number}")
        plt.savefig(fr'venn_plot/{sample}_venn_plot.pdf')
