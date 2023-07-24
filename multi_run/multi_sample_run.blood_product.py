#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
血肿参数，从raw_data开始
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
    parser.add_argument('--nt', '-n', help='qsub投递线程数', type=int)
    return parser.parse_args()


def make_shell(tumor, fq1, fq2, adapter, platform, outdir, sh):
    cmd = rf"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

TUMOR_SM={tumor}
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
PL={platform}

FASTQ_1={fq1}
FASTQ_2={fq2}
ADAPTER={adapter}

FASTA="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/J_hg19_reference/hg19_reference_with_NC.fasta"
KNOWN_DBSNP="/mnt/share02/zhouyj/database/human_b37_reference/dbsnp_138.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NLP.bed"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

NT=48
START_DIR="$PWD/TNscope"
BCFTOOLS_BINARY="/mnt/share01/tools/miniconda/envs/bcftools/bin/bcftools"\

# ******************************************
# Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# Run fastp
# ******************************************
/mnt/share01/tools/bin/fastp \
    -i $FASTQ_1 \
    -I $FASTQ_2 \
    -o $TUMOR_SM.chean.R1.fastq.gz \
    -O $TUMOR_SM.chean.R2.fastq.gz \
    --adapter_fasta $ADAPTER \
    --thread 16 \
    --json $TUMOR_SM.fastp.json \
    --html $TUMOR_SM.fastp.html\
    
python /mnt/share02/xuefei/clinical_WGBS/script/QC/fastp_plot.py \
        --fastp-json $TUMOR_SM.fastp.json \
        --outdir . \
        --samplename $TUMOR_SM

# ******************************************************
# Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************************
$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\\tID:$TUMOR_RGID\\tSM:$TUMOR_SM\\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $TUMOR_SM.chean.R1.fastq.gz $TUMOR_SM.chean.R2.fastq.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $NT --sam2bam -i -
  
# ******************************************
# Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    ${{INTERVAL_FILE:+--interval_padding 50 --interval $INTERVAL_FILE}} \
    -i tumor_sorted.bam \
    --algo TNscope --disable_detector sv --trim_soft_clip \
    --tumor_sample $TUMOR_SM --dbsnp $KNOWN_DBSNP \
    --assemble_mode 4 --disable_detector sv \
    --min_tumor_allele_frac 0.005 --filter_t_alt_frac 0.005 \
    --clip_by_minbq 1 --min_init_tumor_lod 3.0 --min_tumor_lod 3.0 --resample_depth 100000 \
    $TUMOR_SM.output_tnscope.pre_filter.vcf.gz

# ******************************************
# Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate ${{INTERVAL_FILE:+-R $INTERVAL_FILE}} -x "FILTER/triallelic_site" $TUMOR_SM.output_tnscope.pre_filter.vcf.gz | \
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" | \
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \
    $BCFTOOLS_BINARY norm -f $FASTA -m +any | $BCFTOOLS_BINARY sort |
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $TUMOR_SM.output_tnscope.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' $TUMOR_SM.output_tnscope.filtered.vcf.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $TUMOR_SM.output_tnscope.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
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
        for line in f:
            tumor, platform, fq1, fq2 = line.strip().split('\t')
            if platform == 'T7':
                adapter = '/mnt/share01/tools/analysis_module/adapter/T7/adapter.fa'
            else:
                adapter = '/mnt/share01/tools/analysis_module/adapter/illumina/adapter_illumina.fa'
            if not os.path.exists(f'{wd}/{tumor}'):
                os.mkdir(f'{wd}/{tumor}')
            make_shell(tumor, fq1, fq2, adapter, platform, f'{wd}/{tumor}', f'{tumor}.sentieon.sh')
            if argv.qsub:
                do_qsub(f'{wd}/{tumor}/{tumor}.sentieon.sh', argv.nt)
            else:
                print(f'write {wd}/{tumor}/{tumor}.sentieon.sh done.')
