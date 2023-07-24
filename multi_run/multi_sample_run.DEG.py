#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
卡替RNA科研（研发）RNA差异表达项目批量投递diff_gene_Sample_vs_Control.sh（因样本名含有特殊符号重写部分脚本）
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-11-30 15:45:21"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import argparse
import os
import subprocess


def get_args():
    parser = argparse.ArgumentParser(description='批量刷出多个样本的shell并投递运行')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--qsub', '-q', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', '-n', help='qsub投递线程数', type=int, default=1)
    return parser.parse_args()


def make_shell(sample, control, outdir, sh):
    cmd = fr"""#!/bin/sh
set -eo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

source /mnt/share03/tumor_pipeline/Somatic/DNA/miniconda3/bin/activate RNA && \
cat /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/Quantification/{sample}/{sample}.genes.results \
    /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/Quantification/{control}/{control}.genes.results \
| grep -v gene_id | cut -f1 | sort -u > \
    /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/genelist.txt

perl /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/script/extractInfo.pl \
    -i /mnt/share03/tumor_pipeline/Somatic/RNA/RNA_pipeline_v2/database_prepare/hsa/ref_hg19_top_level.format.gff.filter.gtf \
    -g /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/genelist.txt \
    -o /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/geneinfo.txt \
    -l /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/genelength.txt

perl /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/script/EBSeq_v2.pl \
    -i /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/Quantification/merge/merged_readcount \
    -c /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/geneinfo.txt \
    -a R{sample} \
    -b R{control} \
    -n1 R{sample} \
    -n2 R{control} \
    -op /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control \
    -R /mnt/share03/tumor_pipeline/Somatic/DNA/miniconda3/envs/RNA/bin/R && \

perl /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/script/DGEplot_heatmap_volcano_point.pl \
    -p /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control} \
    -c R{sample}:R{control}

rm -rf /mnt/share_fsd/zhouyj/research/RNA/kati_diff_expression_20220913/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/DDN22019843.known_diff_expression_gene_filter.tsv
echo ==== end__ at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)

def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt if nt > 1 else nt + 1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)

    with open(sample_list) as f:
        f.readline()
        for line in f:
            sample, control = line.strip().split('\t')
            make_shell(sample, control, f'{wd}/{sample}_vs_{control}/DiffExpression/Sample_vs_Control', 'diff_gene_Sample_vs_Control.sh')
            if argv.qsub:
                do_qsub(f'{wd}/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/diff_gene_Sample_vs_Control.sh', argv.nt)
            else:
                print(f'write {wd}/{sample}_vs_{control}/DiffExpression/Sample_vs_Control/diff_gene_Sample_vs_Control.sh done.')
