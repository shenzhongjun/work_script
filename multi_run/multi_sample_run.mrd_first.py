#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
Sentieon：mrd单样本bam文件开始，最小阈值call变异
"""

__author__ = "ZhouYiJie"
__email__ = "zhouyijie@zhuanhuayixue.org"
__version__ = "1.0"
__date__ = "2021-11-30 15:45:21"
__copyright__ = "Copyright © 2021 Chigene, All rights reserved"

import os
import subprocess
import argparse


def get_args():
    parser = argparse.ArgumentParser(description='批量刷出多个样本的shell并投递运行')
    parser.add_argument('--sample_list', '-s', help='样本列表', required=True)
    parser.add_argument('--qsub', '-q', help='是否直接投递，默认否', action='store_true')
    parser.add_argument('--nt', '-n', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_shell(tumor, bam, outdir, sh):
    cmd = fr"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
INTERVAL_FILE="/mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/MRD_test/ZRnDZL/mrd.snv.bed"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR="$PWD/TNscope_highcov"
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"
MIN_QUAL=1
MIN_AF=0.00000001
PV=0.99

# ******************************************
# Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    --interval $INTERVAL_FILE --interval_padding 150 \
    -i {bam} \
    --algo TNscope \
    --tumor_sample {tumor} \
    --min_base_qual $MIN_QUAL --clip_by_minbq 1 \
    --min_tumor_allele_frac $MIN_AF --filter_t_alt_frac $MIN_AF \
    --disable_detector sv --trim_soft_clip \
    --max_fisher_pv_active $PV --resample_depth 500000 --assemble_mode 4 \
    {tumor}.force_call.pre_filter.vcf.gz || \
    {{ echo "TNscope failed"; exit 1; }}

# ******************************************
# Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate -R $INTERVAL_FILE -x "FILTER/triallelic_site" {tumor}.force_call.pre_filter.vcf.gz | \
    $BCFTOOLS_BINARY filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" | \
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" | \
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \
    $BCFTOOLS_BINARY sort | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.force_call.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' {tumor}.force_call.filtered.vcf.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - {tumor}.force_call.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt+1} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    tumors = [t.strip() for t in open(sample_list).readlines()]
    bams = [f'{wd}/{t}/bwa_alignment/{t}/{t}.sorted.bam' for t in tumors]
    for t, b in zip(tumors, bams):
        sh = 'run.mrd.sh'
        make_shell(t, b, f'{wd}/{t}', sh)
        do_qsub(f'{wd}/{t}/{sh}', argv.nt) if argv.qsub else print('done.')



