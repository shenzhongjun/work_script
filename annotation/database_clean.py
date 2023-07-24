#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
肿瘤变异注释相关数据库数据清洗
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2023-6-6 16:34:10"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import glob
import os
import time

import pandas as pd


def get_args():
    parser = argparse.ArgumentParser(description='肿瘤变异注释相关数据库数据清洗。')
    parser.add_argument('--hotspot', help='原始癌症热点突变数据库路径')
    parser.add_argument('--osgene', help='原始oncogenes/suppressorgenes原癌基因抑癌基因数据库路径')
    # parser.add_argument('--uniprot', help='原始uniprot蛋白结构域数据库路径')      # 表里坐标是hg38的，改用AA posion匹配
    parser.add_argument('--cbp', help='cBioPortal所有基因注释数据下载结果文件夹')
    parser.add_argument('--ckb', help='CKB所有基因注释数据下载结果文件夹')
    parser.add_argument('--merge', help='合并cBioPortal、CKB和白名单所有位点')
    return parser.parse_args()


def get_output_path(input_path, suffix):
    datenow = time.strftime("%Y%m%d")
    input_path = os.path.abspath(input_path)
    dirname, filename = os.path.split(input_path)
    file, ori_suffix = os.path.splitext(filename)
    return f'{dirname}/{file}.clean.{datenow}.{suffix}' if suffix else f'{dirname}/{file}.{datenow}{ori_suffix}'


def clean_hotspot(hotspot_path):
    """
    把原始数据库snv、indel两个sheet整合到一起，新增AA_Change_Count、Variant_Amino_Acid、AA_Change列
    """
    df = pd.read_excel(hotspot_path, sheet_name=None)  # 读入snv和indel两个sheet，储存在字典df中
    cols = ['Hugo_Symbol', 'Amino_Acid_Position', 'Mutation_Count', 'Reference_Amino_Acid', 'Variant_Amino_Acid',
            'Is_repeat', 'seq', 'length', 'inOncokb', 'Samples']
    snv_df = df['SNV-hotspots'].copy()[cols]
    indel_df = df['INDEL-hotspots'].copy()[cols]
    snv_df['Mutation_type'] = 'snv'
    indel_df['Mutation_type'] = 'indel'
    final_df = snv_df.append(indel_df)
    final_df['AA_Change_Count'] = final_df['Variant_Amino_Acid'].str.split(':', expand=True)[1].astype('int')
    final_df['Variant_Amino_Acid'] = final_df['Variant_Amino_Acid'].str.split(':', expand=True)[0]
    final_df['AA_Change'] = final_df.apply(
        lambda x: f"p.{x['Reference_Amino_Acid'].split(':')[0]}{x['Amino_Acid_Position']}{x['Variant_Amino_Acid']}"
        if x['Mutation_type'] == 'snv' else f"p.{x['Variant_Amino_Acid']}",
        axis=1
    )
    final_df = final_df.reindex(
        ['Hugo_Symbol', 'Amino_Acid_Position', 'Mutation_Count', 'AA_Change', 'AA_Change_Count', 'Mutation_type',
         'Reference_Amino_Acid', 'Is_repeat', 'seq', 'length', 'inOncokb', 'Samples'],
        axis=1
    )
    final_df = final_df.sort_values(['Hugo_Symbol', 'Amino_Acid_Position'])

    print(final_df)
    print(f"Cancer Hotspot数据库输出路径： {get_output_path(hotspot_path, 'txt')}")
    final_df.to_csv(get_output_path(hotspot_path, 'txt'), sep='\t', index=False)


def clean_osgene(osgene_path):
    """原癌基因抑癌基因原始数据处理成方便注释的格式"""
    df = pd.read_excel(osgene_path, sheet_name=3)
    oncogenes = set(df['Oncogene'].dropna().to_list())
    supsgenes = set(df['Suppressor_genes'].dropna().to_list())

    final_df = pd.DataFrame()
    final_df['基因名'] = sorted(list(oncogenes.union(supsgenes)))
    final_df['原癌基因'] = final_df['基因名'].apply(lambda x: '是' if x in oncogenes else '否')
    final_df['抑癌基因'] = final_df['基因名'].apply(lambda x: '是' if x in supsgenes else '否')

    print(final_df)
    print(f"原癌基因抑癌基因数据库输出路径： {get_output_path(osgene_path, 'txt')}")
    final_df.to_csv(get_output_path(osgene_path, 'txt'), sep='\t', index=False)


def clean_uniprot(uniprot_path):
    """已弃用。uniprot数据库转bed"""
    df = pd.read_excel(uniprot_path)
    df = df[['Chromosome', 'Start', 'End', '中文名称']]
    df = df.sort_values(['Chromosome', 'Start'])
    df.to_csv(get_output_path(uniprot_path, 'bed'), sep='\t', index=False, header=False)


def clean_cbioportal(cbp_path):
    """
    cbioportal各基因数据合并后删除冗余
    数据清洗原则：去除fusion数据、去除OncoKB注释为Unknown的无效数据、保留Mut in Sample以做参考（虽然价值不大）
    """
    merged_df = pd.DataFrame()
    read_cols = ['Protein Change', 'Annotation', '# Mut in Sample', 'Mutation Type', 'Chromosome',
                 'Start Pos', 'End Pos', 'Ref', 'Var', 'HGVSg', 'HGVSc', 'COSMIC', 'ClinVar']
    raw_cols = ['Gene', 'Protein Change', 'OncoKB', 'Level', 'Resistance', '# Mut in Sample', 'Mutation Type',
                'Chromosome',
                'Start Pos', 'End Pos', 'Ref', 'Var', 'HGVSg', 'HGVSc', 'COSMIC', 'ClinVar', 'Annotation']
    clean_cols = ['Gene', 'HGVSg', 'Protein Change', 'OncoKB', 'Mutation Type', '# Mut in Sample']
    final_cols = ['Gene', 'HGVSg', 'Protein Change', 'OncoKB', 'Mutation Type', 'Mut in Sample']

    if os.path.exists(fr'{cbp_path}/../redownload_list.txt'):
        os.remove(fr'{cbp_path}/../redownload_list.txt')

    for gene in os.listdir(cbp_path):
        print(gene)
        if not glob.glob(fr'{cbp_path}/{gene}/cBioPortal_*_{gene}.txt'):
            print(f'{gene}文件夹下没有数据！')
            continue
        try:
            gene_df = pd.read_table(glob.glob(fr'{cbp_path}/{gene}/cBioPortal_*_{gene}.txt')[0], dtype='str')
            gene_df = gene_df[read_cols]
            if (gene_df['Annotation'].str.split(';', expand=True)[0].str.split(':', expand=True)[
                    1].str.strip() == 'NA').all():
                print(f'{gene}数据下载不完整，应重新下载！')
                with open(fr'{cbp_path}/../redownload_list.txt', 'a+') as a:
                    a.write(f'{gene}\n')
                continue

            gene_df = gene_df.sort_values(['Chromosome', 'Start Pos']).drop_duplicates()
            gene_df['Gene'] = gene
            tmp = gene_df['Annotation'].str.split(';', expand=True)[0].str.split(':', expand=True)[
                1].str.strip().str.split(',', expand=True)
            gene_df['OncoKB'] = tmp[0]
            gene_df['Level'] = tmp[1].str.split(' ', expand=True)[1]
            gene_df['Resistance'] = tmp[2].str.split(' ', expand=True)[1]

            if merged_df.empty:
                merged_df = gene_df
            else:
                merged_df = pd.concat([merged_df, gene_df], ignore_index=True)
            # if gene == 'MAX': break
        except Exception as ex:
            print(f'{gene}数据处理失败！{ex}')  # 原始数据不含注释列，可忽略

    merged_df = merged_df[raw_cols]
    merged_df.to_csv(fr'{cbp_path}/../cBioPortal.raw.{time.strftime("%Y%m%d")}.txt', sep='\t', index=False)

    # 开始数据清洗
    merged_df = merged_df[clean_cols]
    merged_df = merged_df[(merged_df['Mutation Type'] != 'fusion') & (merged_df['OncoKB'] != 'Unknown')]
    merged_df = merged_df.drop_duplicates()
    merged_df['HGVSg'] = merged_df['HGVSg'].apply(lambda x: 'chr' + x)
    # 处理# Mut in Sample。注：医学最新解释：理解有误，这列不要了
    merged_df = merged_df.rename(columns={'# Mut in Sample': 'Mut in Sample'})
    merged_df['Mut in Sample'] = pd.to_numeric(merged_df['Mut in Sample'])  # 注意：包含空值，无法转为int
    merged_df = merged_df.groupby(['Gene', 'HGVSg', 'Protein Change', 'OncoKB', 'Mutation Type'], as_index=False)[
        'Mut in Sample'].sum()      # 合并相同行，Mut in Sample相加
    # 上一步合并后剪切位点的Protein Change有重复，需单独合并
    splice_df = merged_df[merged_df['Mutation Type'].isin(['Splice_Site', 'Splice_Region'])].copy()
    cds_df = merged_df[~merged_df['Mutation Type'].isin(['Splice_Site', 'Splice_Region'])].copy()
    splice_df = splice_df.groupby(['Gene', 'HGVSg', 'OncoKB', 'Mutation Type'], as_index=False).agg(
        {'Mut in Sample': 'sum', 'Protein Change': lambda x: ';'.join(x)}
    )
    splice_df = splice_df.reindex(final_cols, axis=1)  # 合并后Protein Change变为最后一列，需重新排序
    # 上一步完成后把剪切位点数据和CDS数据重新合并到一起
    merged_df = pd.concat([cds_df, splice_df], ignore_index=True)
    merged_df['Mut in Sample'] = merged_df['Mut in Sample'].astype(int)  # 转为int使输出不带小数点
    # 按Gene名和p.大小排序后删除pos列
    merged_df['pos'] = merged_df['Protein Change'].str.extract(r'.*?(\d+).*').astype(float)     # 有na，转为float
    merged_df = merged_df.sort_values(by=['Gene', 'pos']).drop('pos', axis=1)
    merged_df.to_csv(fr'{cbp_path}/../cBioPortal.clean.{time.strftime("%Y%m%d")}.txt', sep='\t', index=False)
    print(merged_df)
    print(f"cBioPortal数据库输出路径： {cbp_path}/../cBioPortal.clean.{time.strftime('%Y%m%d')}.txt")


def clean_ckb(ckb_path):
    """CKB数据库，可直接合并"""
    merged_df = pd.DataFrame()

    for gene_file in os.listdir(ckb_path):
        gene_name = os.path.splitext(gene_file)[0]
        print(gene_name)
        gene_df = pd.read_excel(fr'{ckb_path}/{gene_file}', dtype='str')
        gene_df['Gene'] = gene_name
        if merged_df.empty:
            merged_df = gene_df
        else:
            merged_df = pd.concat([merged_df, gene_df], ignore_index=True)
    merged_df.to_csv(fr'{ckb_path}/../ckb.raw.{time.strftime("%Y%m%d")}.txt', sep='\t', index=False)
    # data clean
    merged_df = merged_df[['Gene', 'Alteration', 'protein effect']]
    merged_df = merged_df[merged_df['protein effect'] != 'unknown']
    merged_df = merged_df[merged_df['Alteration'].str.contains(r'\d')]     # 去除cnv信息
    merged_df.to_csv(fr'{ckb_path}/../ckb.clean.{time.strftime("%Y%m%d")}.txt', sep='\t', index=False)
    print(merged_df)
    print(f'CKB数据库输出路径： {ckb_path}/../ckb.clean.{time.strftime("%Y%m%d")}.txt')


def merge_dbs(dbs):
    cbp, ckb, white = dbs.split(',')
    cbp_df = pd.read_table(cbp, dtype=str)
    cbp_df['merge_id'] = cbp_df['Gene'].str.cat(cbp_df['Protein Change'], sep='_')

    ckb_df = pd.read_table(ckb)
    ckb_df = ckb_df.rename(columns={'Gene': 'ckbGene'})
    ckb_df['merge_id'] = ckb_df['ckbGene'].str.cat(ckb_df['Alteration'], sep='_')

    white_df = pd.read_excel(white)
    white_df = white_df[white_df['变异类型(除“特殊检测”外其他变异均精确匹配氨基酸，特殊检测需识别给定区域范围内的特定变异类型)'] != '特殊检测']
    white_df = white_df[['基因', '变异', '临床意义']]
    white_df = white_df.dropna(subset=['临床意义'], axis=0)
    white_df['merge_id'] = white_df['基因'].str.cat(white_df['变异'], sep='_')

    merge_df = cbp_df.merge(white_df, on='merge_id', how='outer')
    merge_df['Gene'] = merge_df['Gene'].fillna(merge_df['基因'])
    merge_df['Protein Change'] = merge_df['Protein Change'].fillna(merge_df['变异'])
    merge_df['OncoKB'] = merge_df['OncoKB'].fillna(merge_df['临床意义'])
    merge_df = merge_df.drop(columns=['基因', '变异', '临床意义'])

    merge_df = merge_df.merge(ckb_df, on='merge_id', how='outer')
    merge_df['Gene'] = merge_df['Gene'].fillna(merge_df['ckbGene'])
    merge_df['Protein Change'] = merge_df['Protein Change'].fillna(merge_df['Alteration'])
    merge_df['OncoKB'] = merge_df['OncoKB'].fillna(merge_df['protein effect'])
    merge_df = merge_df.drop(columns=['ckbGene', 'Alteration', 'protein effect', 'merge_id'])

    # 按Gene名和p.大小排序后删除pos列
    merge_df['pos'] = merge_df['Protein Change'].str.extract(r'.*?(\d+).*').astype(float)  # 有na，转为float
    merge_df = merge_df.sort_values(by=['Gene', 'pos']).drop('pos', axis=1).fillna('-')

    merge_df.to_csv(fr'{os.path.dirname(cbp)}/../cbp_ckb_white.merge.{time.strftime("%Y%m%d")}.txt', sep='\t', index=False)
    print(merge_df)
    print(f'cBioPortal、CKB、白名单列表合并输出路径： {os.path.dirname(cbp)}/../cbp_ckb_white.merge.{time.strftime("%Y%m%d")}.txt')


if __name__ == "__main__":
    args = get_args()
    if args.hotspot:
        clean_hotspot(args.hotspot)
    if args.osgene:
        clean_osgene(args.osgene)
    # if args.uniprot:
    #     clean_uniprot(args.uniprot)
    if args.cbp:
        clean_cbioportal(args.cbp)
    if args.ckb:
        clean_ckb(args.ckb)
    if args.merge:
        merge_dbs(args.merge)
