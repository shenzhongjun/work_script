#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
MRD位点筛选脚本，运行pyclone后依据规则挑选指定个数的MRD位点。
待更新：
加入靶药位点筛选；与pyclone脚本合并；结果格式优化，删除merge列和pyclone列
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "0.1.2"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import pandas as pd
import subprocess


def get_args():
    parser = argparse.ArgumentParser(description='MRD位点筛选脚本，从merge表出发，依据规则挑选指定个数的MRD位点。')
    parser.add_argument('--merge', help='merge表路径', required=True)
    parser.add_argument('--cnvkit', help="cnvkit结果路径", required=True)
    parser.add_argument('--purity', help='肿瘤纯度', required=True)
    parser.add_argument('--target_results', help='检测结果汇总路径', required=True)
    parser.add_argument('--wes_tumor_id', help='WES样本名', required=True)
    parser.add_argument('--clonal', help='克隆性造血数据库路径',
                        default='/home/zhouyj/script/mrd/克隆性造血突变库-20220421.xlsx')
    parser.add_argument('--wes_stats', help='NT01T芯片突变频率统计数据库',
                        default='/home/zhouyj/script/mrd/NT01T.somatic_stats.xls')
    parser.add_argument('--site_num', help='要筛选的位点数目', default=50)
    parser.add_argument('--order_id', help='订单编号', required=True)
    parser.add_argument('--out', help='结果输出目录', required=True)
    return parser.parse_args()


def get_cnv(cnvkit_df, x):
    """
    从cnvkit结果中获取突变位点的cnv值
    chromosome      start   end     gene    log2    cn      depth   p_ttest probes  weight
    """
    pd_cnvkit_fil = cnvkit_df[
        (cnvkit_df['chromosome'] == x['染色体名称']) &
        (cnvkit_df['start'] <= x['起始位置']) &
        (cnvkit_df['end'] >= int(x['终止位置']))]
    pd_cnvkit_fil.index = range(len(pd_cnvkit_fil))
    if pd_cnvkit_fil.empty:
        return '0'
    elif len(pd_cnvkit_fil) == 1:
        return pd_cnvkit_fil.loc[0, 'cn']
    else:
        print('sth is wrong with cnv number')  # f"{chrom}:{start}-{end}"
        return '0'


def run_pyclone(pd_merge, wes_id, purity, outdir):
    """
    获得最终的Pyclone-vi 的input
    mutation_id	sample_id	ref_counts	alt_counts	normal_cn	major_cn	minor_cn	tumour_content
    CRUK0001:1:1564541:C	R1	178	0	2	2	1	0.21
    """
    out_columns = ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts',
                   'normal_cn', 'major_cn', 'minor_cn', 'tumour_content']
    pd_merge['sample_id'] = wes_id
    pd_merge['normal_cn'] = 2
    pd_merge['minor_cn'] = 0
    pd_merge['tumour_content'] = purity
    pd_out = pd_merge[['突变编号', 'sample_id', f"总深度(待定,{wes_id})", f"突变深度(待定,{wes_id})",
                      'normal_cn', 'major_cn', 'minor_cn', 'tumour_content']]
    pd_out.columns = out_columns
    pd_out = pd_out.drop_duplicates(subset=['mutation_id'], keep='first')
    pd_out.to_csv(f"{outdir}/{wes_id}.pyclone-vi.txt", sep="\t", index=False)
    cmd = fr"""
/mnt/share02/lixx/02.software/envs_python3.7/PyClone-VI/bin/pyclone-vi fit \
    -i {wes_id}.pyclone-vi.txt \
    -o {wes_id}.pyclone-vi.h5 \
    -c 40 -d beta-binomial -r 10

/mnt/share02/lixx/02.software/envs_python3.7/PyClone-VI/bin/pyclone-vi write-results-file \
    -i {wes_id}.pyclone-vi.h5 \
    -o {wes_id}.pyclone-vi.output.txt
"""
    subprocess.call(cmd, shell=True)
    pyclone_df = pd.read_table(f"{outdir}/{wes_id}.pyclone-vi.output.txt", low_memory=False)
    pyclone_df = pyclone_df.sort_values('cellular_prevalence', ascending=False)     # 排序后取细胞占比大的id作为主克隆id
    main_cluster = pyclone_df.loc[0, 'cluster_id']
    pyclone_df = pyclone_df[pyclone_df['cluster_id'] == main_cluster]
    return pyclone_df


if __name__ == "__main__":
    args = get_args()
    wes_id = args.wes_tumor_id
    site_num = int(args.site_num)
    out_dir = args.out
    conclusion = f'结论(待定,{wes_id})' if wes_id else '结论'
    merge_df = pd.read_table(args.merge, low_memory=False, dtype={'终止位置': 'str'})
    clonal_df = pd.read_excel(args.clonal, dtype={'Stop': 'str'})
    wes_df = pd.read_table(args.wes_stats, low_memory=False)
    wes_df = wes_df[(wes_df['频率'] >= 0.05)]   # 获取WES高频变异，阈值为5%
    cnvkit_df = pd.read_csv(args.cnvkit, sep='\t')

    merge_df['merge'] = merge_df['染色体名称'].str.cat(merge_df['终止位置'], sep='_')
    clonal_df['merge'] = clonal_df['Chr'].str.cat(clonal_df['Stop'], sep='_')

    merge_df.dropna(axis=0, how='any', subset=[f"突变深度(待定,{wes_id})"])     # 去空值
    merge_df.drop_duplicates(subset='突变编号', keep='first', ignore_index=True)  # 必须去重否则pyclone报错

    merge_df = merge_df[
        (merge_df[f"突变深度(待定,{wes_id})"] >= 8) &
        (merge_df['突变可靠性'] == '相对可靠') &
        (~merge_df[conclusion].isin(['遗传性突变', '突变来源不确定'])) &
        (~merge_df['突变类型'].isin(['同义突变', '非编码区突变', '剪切位点附近的突变'])) &
        ((merge_df['1000genomeMAF'].isna()) | (merge_df['1000genomeMAF'] <= 0.01)) &
        ((merge_df['ExAC-eas'].isna()) | (merge_df['ExAC-eas'] <= 0.01)) &
        ((merge_df['ExAC-biggest'].isna()) | (merge_df['ExAC-biggest'] <= 0.01)) &
        ((merge_df['ESP_MAF'].isna()) | (merge_df['ESP_MAF'] <= 0.01)) &
        ((merge_df['gnomad_eas_maf'].isna()) | (merge_df['gnomad_eas_maf'] <= 0.01))]
    print(f'step1:筛选merge表\n剩余变异{merge_df.shape[0]}个\n')

    merge_df = merge_df[~merge_df['突变编号'].isin(wes_df['突变编号'].tolist())]
    print(f'step2:过滤WES芯片高频变异\n剩余变异{merge_df.shape[0]}个\n')

    merge_df = merge_df[~merge_df['merge'].isin(clonal_df['merge'].tolist())]
    print(f'step3:过滤克隆性造血数据库\n剩余变异{merge_df.shape[0]}个\n')

    merge_df['major_cn'] = merge_df.apply(lambda x: get_cnv(cnvkit_df, x), axis=1)
    pyclone_df = run_pyclone(merge_df, wes_id, args.purity, out_dir)

    merge_df2 = merge_df.copy()     # 筛选主克隆突变和靶点库突变前的所有变异
    merge_df = merge_df[merge_df['突变编号'].isin(pyclone_df['mutation_id'].tolist())]
    merge_df['MRD来源'] = '主克隆突变'
    print(f"step4:筛选主克隆突变\n剩余变异{merge_df.shape[0]}个\n")

    if merge_df.shape[0] < site_num:
        targeted_df = pd.read_table(args.target_results)
        targeted_df = targeted_df[targeted_df['类型'] == '靶向药物相关突变']
        i = 0
        for gene, pchange in zip(targeted_df['变异基因'].tolist(),
                                 targeted_df['检测结果'].str.split(' ', expand=True)[1].tolist()):
            tmp = merge_df2[(merge_df2['基因'] == gene) & (merge_df2['氨基酸改变'] == pchange)]
            if tmp['突变编号'].iloc[0] in merge_df['突变编号'].tolist():
                print(f'{gene}:{pchange} 已经包含在主克隆突变内，跳过\n')
            else:
                merge_df = merge_df.append(tmp)
                i += 1
        print(f'step5:突变个数不足{site_num}个，加入靶点变异{i}个\n')

        if merge_df.shape[0] < site_num:
            need_count = site_num - merge_df.shape[0]
            merge_df2 = merge_df2[(~merge_df2['突变编号'].isin(merge_df['突变编号'].tolist())) &
                                  (merge_df2[f'突变深度(待定,{wes_id})'] >= 8)]
            merge_df2 = merge_df2.sort_values([f'突变率(待定,{wes_id})'], ascending=False)
            merge_df2 = merge_df2.iloc[0:need_count]
            merge_df2['MRD来源'] = '高频体细胞突变'
            merge_df = merge_df.append(merge_df2)
            print(f'step6:突变个数依然不足{site_num}个，加入高频体细胞变异{need_count}个\n')

        merge_df = merge_df.sort_values(['MRD来源', f'突变率(待定,{wes_id})'], ascending=[True, False])
    else:
        merge_df = merge_df.sort_values(f'突变率(待定,{wes_id})', ascending=False)
        merge_df = merge_df.iloc[0:50]
        print(f"主克隆突变大于等于{site_num}个，按突变率从高到低取前50个\n")

    merge_df = merge_df.drop(['merge', 'normal_cn', 'minor_cn', 'tumour_content', 'major_cn', 'sample_id'], axis=1)
    merge_df.to_csv(f'{out_dir}/{args.order_id}_mrd_sites.xls', sep='\t', index=False)

