#!/mnt/share01/tools/bin/python
# -*- coding: UTF-8 -*-

import argparse
import os

import pandas as pd


def Init_args():
    parser = argparse.ArgumentParser(description="整合cnvkit和merge表生成Pyclone-VI的input结果")
    parser.add_argument('--merge', help="the merge file from pipline", required=True)
    parser.add_argument('--cnvkit', help="the result generated by cnvkit", required=True)
    parser.add_argument('--tissue_id', help='the sample id of tissue', required=True)
    parser.add_argument('--purity', help='the purity of tissue', required=True)
    parser.add_argument('--outdir', help='the output directory', default=os.getcwd(), required=False)
    parser.add_argument('--wes_stats', help='NT01T芯片突变频率统计数据库', default='/home/zhouyj/script/mrd/NT01T.somatic_stats.xls')
    args = parser.parse_args()
    return args


def Merge_filter(merge_file, tissue_id, wes_stats):
    """
    merge表格进行过滤
    """
    wes_df = pd.read_table(wes_stats, low_memory=False)
    wes_df = wes_df[wes_df['频率'] >= 0.05]   # 获取WES高频变异，阈值为5%
    pd_merge = pd.read_csv(merge_file, sep='\t', low_memory=False)
    pd_merge = pd_merge.dropna(axis=0, how='any', subset=[f"突变深度(待定,{tissue_id})"])     # 去空值
    pd_merge = pd_merge.drop_duplicates(subset='突变编号', keep='first', ignore_index=True)  # 必须去重否则pyclone报错
    pd_merge = pd_merge[pd_merge[f"突变深度(待定,{tissue_id})"] >= 8]
    pd_merge = pd_merge[
        (~pd_merge[f"结论(待定,{tissue_id})"].isin(['遗传性突变', '突变来源不确定'])) &     # , '突变来源不确定'
        (pd_merge['突变可靠性'] == '相对可靠') &
        (~pd_merge['突变类型'].isin(['同义突变', '非编码区突变', '剪切位点附近的突变']))]
    pd_merge = pd_merge[~pd_merge['突变编号'].isin(wes_df['突变编号'].tolist())]
    print(f'过滤WES高频突变{pd_merge[pd_merge["突变编号"].isin(wes_df["突变编号"].tolist())].shape[0]}个')
    return pd_merge


def Get_cnv(cnvkit, x):
    """
    从cnvkit结果中获取突变位点的cnv值
    chromosome      start   end     gene    log2    cn      depth   p_ttest probes  weight
    """
    pd_cnvkit = pd.read_csv(cnvkit, sep='\t')
    pd_cnvkit_fil = pd_cnvkit[
        (pd_cnvkit['chromosome'] == x['染色体名称']) &
        (pd_cnvkit['start'] <= x['起始位置']) &
        (pd_cnvkit['end'] >= x['终止位置'])]
    pd_cnvkit_fil.index = range(len(pd_cnvkit_fil))
    if pd_cnvkit_fil.empty:
        return '0'
    elif len(pd_cnvkit_fil) == 1:
        return pd_cnvkit_fil.loc[0, 'cn']
    else:
        print('sth is wrong with cnv number')  # f"{chrom}:{start}-{end}"
        return '0'


def Get_out(pd_merge, tissue_id, purity, outdir):
    """
    获得最终的Pyclone-vi 的input
    mutation_id	sample_id	ref_counts	alt_counts	normal_cn	major_cn	minor_cn	tumour_content
    CRUK0001:1:1564541:C	R1	178	0	2	2	1	0.21
    """
    out_columns = ['mutation_id', 'sample_id', 'ref_counts', 'alt_counts', 'normal_cn', 'major_cn', 'minor_cn',
                   'tumour_content']
    pd_merge['sample_id'] = tissue_id
    pd_merge['normal_cn'] = 2
    pd_merge['minor_cn'] = 0
    pd_merge['tumour_content'] = purity
    out_title = ['突变编号', 'sample_id', f"总深度(待定,{tissue_id})", f"突变深度(待定,{tissue_id})", 'normal_cn', 'major_cn',
                 'minor_cn', 'tumour_content']
    pd_out = pd_merge[out_title]
    pd_out.columns = out_columns
    pd_out = pd_out.drop_duplicates(subset=['mutation_id'], keep='first')
    pd_out.to_csv(f"{outdir}/{tissue_id}.pyclone-vi.txt", sep="\t", index=False)


def main():
    argv = Init_args()
    pd_merge = Merge_filter(argv.merge, argv.tissue_id, argv.wes_stats)
    pd_merge['major_cn'] = pd_merge.apply(lambda x: Get_cnv(argv.cnvkit, x), axis=1)
    Get_out(pd_merge, argv.tissue_id, argv.purity, argv.outdir)


if __name__ == "__main__":
    main()