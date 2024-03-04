#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CAP补充实验，从结果文件夹中提取需要的结果
TODO qc合并、wes tmb计算，msi提取
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__date__ = "2024-1-30 14:15:01"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import argparse
import pandas as pd

pd.set_option('display.max_columns', None)

def get_args():
    parser = argparse.ArgumentParser(description='CAP补充实验，从结果文件夹中提取需要的结果')
    parser.add_argument('--sample_list', '-s', help='订单、样本、突变对应信息', required=True)
    parser.add_argument('--dirs', '-d', help='要搜索的文件夹列表,用空格分隔', required=True)
    parser.add_argument('--chip', '-c', help='芯片类型', choices=['npc69', 'wes'], required=True)
    parser.add_argument('--out', '-o', help='输出文件路径', required=True)
    # parser.add_argument('--qc', '-q', help='质控结果输出文件夹', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    argv = get_args()
    df = pd.read_table(argv.sample_list)
    if argv.chip == 'npc69':
        df = df[df['Chip'] == '69panel']
    else:
        df = df[df['Chip'] == 'WES']
    outfile = open(argv.out, 'w')

    for i in argv.dirs.strip().split(' '):
        for index, row in df.iterrows():
            # qc结果使用Linux命令手动合并更便捷
            # if argv.chip == 'npc69':
            #     path1 = os.path.realpath(
            #         f'{i}/{row["Library"]}/anno/{row["Library"]}/01.annotools/{row["Library"]}_NCPLt-NPC69_rmdup_QC.xlsx')
            #     path2 = f'{argv.qc}/{row["Library"]}_NCPLt-NPC69_rmdup_QC.xlsx'
            # else:
            #     path1 = os.path.realpath(
            #         f'{i}/{row["Library"]}/anno/{row["Library"]}/01.annotools/{row["Library"]}_NCPLt-NPC69_rmdup_QC.xlsx')
            #     path2 = f'{argv.qc}/{row["Library"]}_NCPLt-NPC69_rmdup_QC.xlsx'
            #
            # if os.path.exists(path1):
            #     os.system(f'ln -sf {path1} {path2}')
            # else:
            #     if os.path.exists(f'{i}/{row["Library"]}'):
            #         # print(f'{i}/{row["Library"]}没有找到QC文件，请手动核查！')
            #         pass

            if os.path.exists(f'{i}/{row["Library"]}/single_panel/new/{row["Order"]}_merge.txt'):
                merge_new_file = f'{i}/{row["Library"]}/single_panel/new/{row["Order"]}_merge.txt'
            elif os.path.exists(f'{i}/{row["Library"]}/lung/new/{row["Order"]}_merge.txt'):
                merge_new_file = f'{i}/{row["Library"]}/lung/new/{row["Order"]}_merge.txt'
            elif os.path.exists(f'{i}/{row["Library"]}/single_wes/mutation.txt.civic.protein_merge'):
                merge_new_file = f'{i}/{row["Library"]}/single_wes/mutation.txt.civic.protein_merge'
                print(merge_new_file)
            else:
                merge_new_file = ''
                continue
            if pd.isnull(row['OriginalSomatic']): continue
            print(row['OriginalSomatic'].split(';'))
            origi_muts = ['@'.join(i.split(' ')[0:3]) for i in row['OriginalSomatic'].split(';')]
            # 读取merge表
            merge_df = pd.read_table(merge_new_file)
            merge_df = merge_df.rename(columns={merge_df.columns[3]: '突变率'})
            # 在此处计算TMB
            tmb_mutation_frame = merge_df.copy()
            tmb_mutation_frame = tmb_mutation_frame[
                (tmb_mutation_frame["突变可靠性"] == "相对可靠") &
                (tmb_mutation_frame["结论"].isin(["体细胞较可靠突变,可能无害", "体细胞较可靠突变,可能有害"])) &
                (~tmb_mutation_frame["突变类型"].isin(["非编码区突变", "剪切位点附近的突变", "剪切位点突变", "同义突变"])) &
                (tmb_mutation_frame["突变率"] >= 0.05)]

            tmb_mutation_frame = tmb_mutation_frame[
                # (merge_df[f"突变深度(待定,{sample_id})"] >= 8) &
                ((tmb_mutation_frame['1000genomeMAF'].isna()) | (tmb_mutation_frame['1000genomeMAF'] <= 0.01)) &
                ((tmb_mutation_frame['ExAC-eas'].isna()) | (tmb_mutation_frame['ExAC-eas'] <= 0.01)) &
                ((tmb_mutation_frame['ExAC-biggest'].isna()) | (tmb_mutation_frame['ExAC-biggest'] <= 0.01)) &
                ((tmb_mutation_frame['ESP_MAF'].isna()) | (tmb_mutation_frame['ESP_MAF'] <= 0.01)) &
                ((tmb_mutation_frame['gnomad_eas_maf'].isna()) | (tmb_mutation_frame['gnomad_eas_maf'] <= 0.01))
                ]
            tmb_mutation_frame['疾病OMIM'] = ''
            tmb_mutation_frame['疾病表型名称'] = ''
            tmb_mutation_frame['遗传方式'] = ''
            tmb_mutation_frame = tmb_mutation_frame.drop_duplicates()
            tmb_num = tmb_mutation_frame.shape[0]
            tmb_value = tmb_num / 33
            tmb_mutation_frame.to_csv(f'tmp/tmb_mutation_{row["Library"]}.txt', sep='\t', index=False)
            print(f'{tmb_value:.2f}')
            df.loc[index, 'TMB'] = f'{tmb_value}'

            merge_df = merge_df.iloc[:, [3, 9, 17, 18]]
            merge_df.drop_duplicates(keep='first', inplace=True)

            if argv.chip == 'npc69':
                merge_df = merge_df[merge_df['突变率'] > 0.005]
            else:
                merge_df = merge_df[merge_df['突变率'] > 0.01]
            merge_df['突变率'] = merge_df['突变率'].apply(lambda x: f'{float(x):.2%}')
            merge_df['外显子'] = merge_df['核酸改变'].str.extract(r'\((.*?)\)')
            merge_df['核苷酸'] = merge_df['核酸改变'].str.replace(r'\(exon\d+\)', '')
            merge_df.loc[merge_df['外显子'].str.startswith('IVS'), '氨基酸改变'] = merge_df['核苷酸']
            merge_df['check'] = merge_df['基因'] + '@' + merge_df['外显子'] + '@' + merge_df['氨基酸改变']
            merge_df['匹配'] = merge_df[merge_df['check'].isin(origi_muts)][['基因', '外显子', '氨基酸改变', '突变率']].apply(lambda row: '@'.join(row), axis=1)
            merge_df = merge_df[merge_df['匹配'].notna()]
            # print(merge_df)
            match_list = merge_df['匹配'].tolist()
            result_string = ';'.join([' '.join(element.split('@')) for element in match_list])
            df.loc[index, 'Somatic'] = result_string
            print(result_string)

            # 在此处抓取MSI结果
            if os.path.exists(f'{i}/{row["Library"]}/msi/{row["Library"]}/{row["Library"]}.msi_summary.txt'):
                msi_file = f'{i}/{row["Library"]}/msi/{row["Library"]}/{row["Library"]}.msi_summary.txt'
                with open(msi_file) as f:
                    f.readline()
                    df.loc[index, 'MSI'] = f.readline().strip().split('\t')[1]
                    print(df.loc[index, 'MSI'])
            elif argv.chip != 'npc69':
                print(f'{row["Library"]}的MSI运行失败，请检查！')
            # break
    print(df)
    df.to_csv(outfile, sep='\t', index=False)









