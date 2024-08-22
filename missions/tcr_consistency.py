#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
生物学重复间TCR克隆一致性统计
"""
__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-4-1 15:13:45"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations

pd.options.display.max_columns = 16


# 计算相似性
def calculate_similarity(df1, df2):
    # 提取两个 DataFrame 中共同的 barcode 数据
    common_barcodes = set(df1['barcode']).intersection(df2['barcode'])
    intersection_count = len(df1.merge(df2, on=['barcode', 'aaSeqCDR3_TRA', 'aaSeqCDR3_TRB']))
    # print(df1.merge(df2, on=['barcode', 'aaSeqCDR3_TRA', 'aaSeqCDR3_TRB']))
    # 计算相似性
    similarity = intersection_count / len(common_barcodes) if len(common_barcodes) != 0 else 0
    return similarity


for number in ['1', '2', '3', '4']:
    dfs = []
    keys = []
    for letter in ['A', 'B', 'C', 'D', 'E', 'F', 'G']:
        df = pd.read_table(glob.glob(fr'{letter}{number}/05.count_vdj/outs/{letter}{number}_*.merge_vdj_clonetype_scRNACluster.xls')[0], dtype='str')
        df = df.loc[
            (~df['aaSeqCDR3_TRA'].isna()) & (~df['aaSeqCDR3_TRB'].isna()) & (df['celltype'] == 'T cells')]
        df['sample'] = f'{letter}{number}'
        dfs.append(df[['sample', 'barcode', 'aaSeqCDR3_TRA', 'aaSeqCDR3_TRB']])
        keys.append(f'{letter}{number}')
    # 比较两两 DataFrame 之间的相似性
    similarity_df = pd.DataFrame(index=keys, columns=keys)
    for df1, df2 in combinations(dfs, 2):
        similarity = calculate_similarity(df1, df2)
        sample1 = df1['sample'].tolist()[0]
        sample2 = df2['sample'].tolist()[0]
        similarity_df.loc[sample1, sample1] = 1
        similarity_df.loc[sample2, sample2] = 1
        similarity_df.loc[sample1, sample2] = similarity
        similarity_df.loc[sample2, sample1] = similarity
    print(similarity_df)
    similarity_df.to_csv(f'similarity/sample{number}_similarity_matrix.txt', sep='\t', index=True, header=True)

    # 使用 Seaborn 绘制相似性热图
    similarity_df = similarity_df.astype(float)
    plt.figure(figsize=(10, 8))
    sns.heatmap(similarity_df, annot=True, fmt='.2f', cmap='YlGnBu')   # , mask=mask, coolwarm
    plt.xlabel('Samples')
    plt.ylabel('Samples')
    plt.title('Similarity Heatmap')
    plt.savefig(f'similarity/sample{number}_similarity_heatmap.png')
    plt.show()

