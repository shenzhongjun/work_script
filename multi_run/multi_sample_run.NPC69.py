#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
69panel从bam开始运行tumor_only测试
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2022-2-24 15:03:52"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import subprocess
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='69panel从bam开始运行tumor_only测试')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--qsub', '-q', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', '-t', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_panel_t_only_shell(tumor, t_bam, outdir, panel, sh):
    bed = '/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC69.bed' if panel == '69' \
        else '/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC980.bed'
    cmd = fr"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

TUMOR_SM={tumor}
TUMOR_BAM={t_bam}

FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
KNOWN_DBSNP="/home/zhouyj/development/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"
INTERVAL_FILE={bed}

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR={outdir}
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# 5. Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    ${{INTERVAL_FILE:+--interval_padding 50 --interval $INTERVAL_FILE}} \
    -i $TUMOR_BAM \
    --algo TNscope --disable_detector sv \
    --tumor_sample $TUMOR_SM --dbsnp $KNOWN_DBSNP \
    --assemble_mode 4 --disable_detector sv --sv_mask_ext 10 \
    --min_tumor_allele_frac 0.008 --filter_t_alt_frac 0.009 \
    --clip_by_minbq 1 --min_init_tumor_lod 3.0 --min_tumor_lod 3.0 \
    --assemble_mode 4 --resample_depth 100000 \
    {tumor}.tumor_only.pre_filter.vcf.gz

# ******************************************
# 6. Variant filtration
# ******************************************
$BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" {tumor}.tumor_only.pre_filter.vcf.gz | \
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \
    $BCFTOOLS_BINARY sort | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.tumor_only.filtered.vcf.gz

perl /mnt/share01/tools/annovar/table_annovar.pl \
    {tumor}.tumor_only.filtered.vcf.gz \
    /mnt/share01/tools/annovar/humandb/ -buildver hg19 \
    -out {tumor}.tumor_only.filtered.anno \
    -remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all \
    -operation g,r,f,f,f,f,f,f -vcfinput

# ******************************************
# 7. Germline Variant filtration
# ******************************************\
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    -i $TUMOR_BAM \
    --algo Haplotyper \
    -d $KNOWN_DBSNP \
    {tumor}.output_hc.vcf.gz

perl /mnt/share01/tools/annovar/table_annovar.pl \
    {tumor}.output_hc.vcf.gz \
    /mnt/share01/tools/annovar/humandb/ -buildver hg19 \
    -out {tumor}.output_hc.anno  \
    -remove -protocol refGene,cytoBand,clinvar_20200316,cosmic68,1000g2015aug_all,snp138,esp6500si_all,ljb26_all \
    -operation g,r,f,f,f,f,f,f -vcfinput

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt + 1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    samples = []
    with open(sample_list) as f:
        f.readline()
        for line in f:
            samples.append(line.strip().split('\t'))
    for i in samples:
        order, tumor, panel = i
        make_panel_t_only_shell(tumor,
                                f'{wd}/{order}_{panel}/{tumor}/{tumor}.final.bam',
                                f'{wd}/{order}_{panel}/{tumor}',
                                f'{panel}',
                                f'run.NPC{panel}.tumor_only.sh')
        if argv.qsub:
            do_qsub(f'{wd}/{order}_{panel}/{tumor}/run.NPC{panel}.tumor_only.sh', argv.nt)
        else:
            print(f'{wd}/{order}_{panel}/{tumor}/run.NPC{panel}.tumor_only.sh done.')
