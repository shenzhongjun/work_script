#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
Sentieon：mrd单样本raw_data文件开始，加入given参数call变异
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
    parser.add_argument('--nt', '-n', help='qsub投递线程数', type=int, default=48)
    return parser.parse_args()


def make_shell(tumor, fq1, fq2, outdir, sh, nt):
    cmd = fr"""#!/bin/sh
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

TUMOR_SM={tumor}
TUMOR_RGID="rg_$TUMOR_SM"
PL="DNBSEQ"
TUMOR_FASTQ_1={fq1}
TUMOR_FASTQ_2={fq2}

FASTA_DIR="/mnt/share02/zhouyj/database/human_b37_reference"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.chr.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.chr.vcf.gz"
KNOWN_INDELS="$FASTA_DIR/1000G_phase1.indels.b37.chr.vcf.gz"
KNOWN_MILLS="$FASTA_DIR/Mills_and_1000G_gold_standard.indels.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/MRD/multiPCR/mrd.snv.bed"
TRUTH_VCF="/mnt/share05/clinical_project/projects/blood_tumor/test/zhouyj_test/MRD/lod/mrd.truth.vcf.gz"

SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

# NT=$(nproc)
NT={nt}
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
# Adapter procsessing
# ******************************************
/mnt/share01/tools/bin/fastp \
    -i $TUMOR_FASTQ_1 \
    -I $TUMOR_FASTQ_2 \
    -o $TUMOR_SM.chean.R1.fastq.gz \
    -O $TUMOR_SM.chean.R2.fastq.gz \
    --adapter_fasta /mnt/share01/tools/analysis_module/adapter/T7/adapter.fa \
    --thread $NT \
    --json $TUMOR_SM.fastp.json \
    --html $TUMOR_SM.fastp.html \

python /mnt/share02/xuefei/clinical_WGBS/script/QC/fastp_plot.py \
    --fastp-json $TUMOR_SM.fastp.json \
    --outdir . \
    --samplename $TUMOR_SM

# ******************************************
# Mapping reads with BWA-MEM, sorting for tumor sample
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" \
    -t $NT -K 10000000 $FASTA $TUMOR_SM.chean.R1.fastq.gz $TUMOR_SM.chean.R2.fastq.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -o tumor_sorted.bam -t $NT --sam2bam -i -
  
# ******************************************
# Somatic and Structural variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT \
    --interval $INTERVAL_FILE --interval_padding 150 \
    -i tumor_sorted.bam \
    --algo TNscope \
    --tumor_sample $TUMOR_SM \
    --min_base_qual $MIN_QUAL --clip_by_minbq 1 \
    --min_tumor_allele_frac $MIN_AF --filter_t_alt_frac $MIN_AF \
    --disable_detector sv --trim_soft_clip \
    --max_fisher_pv_active $PV --resample_depth 500000 --assemble_mode 4 \
    --given $TRUTH_VCF \
    $TUMOR_SM.force_call.pre_filter.vcf.gz || \
    {{ echo "TNscope failed"; exit 1; }}

# ******************************************
# Variant filtration
# ******************************************
$BCFTOOLS_BINARY annotate -R $INTERVAL_FILE -x "FILTER/triallelic_site" $TUMOR_SM.force_call.pre_filter.vcf.gz | \
    $BCFTOOLS_BINARY filter -m + -s "insignificant" -e "(PV>0.25 && PV2>0.25) || (INFO/STR == 1 && PV>0.05)" | \
    $BCFTOOLS_BINARY filter -m + -s "low_qual" -e "QUAL < 10" | \
    $BCFTOOLS_BINARY filter -m + -s "short_tandem_repeat" -e "RPA[0]>=10" | \
    $BCFTOOLS_BINARY filter -m + -s "read_pos_bias" -e "FMT/ReadPosRankSumPS[0] < -5" | \
    $BCFTOOLS_BINARY filter -m + -s "low_af" -e "(FMT/AD[0:1]) / (FMT/AD[0:0] + FMT/AD[0:1]) < 0.0003" | \
    $BCFTOOLS_BINARY sort | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $TUMOR_SM.force_call.filtered.vcf.gz

$BCFTOOLS_BINARY filter -i 'FILTER=="PASS"' $TUMOR_SM.force_call.filtered.vcf.gz | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util vcfconvert - $TUMOR_SM.force_call.filtered.PASS.vcf.gz

echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt } {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    sample_dic = {}
    with open(sample_list) as f:
        for line in f:
            tumor, fq1, fq2 = line.strip().split('\t')
            sample_dic[tumor] = {'tumor': tumor, 'fq1': fq1, 'fq2': fq2}

    for i in sample_dic:
        i = sample_dic[i]
        if not os.path.exists(f"{wd}/{i['tumor']}"):
            os.mkdir(f"{wd}/{i['tumor']}")
        make_shell(i['tumor'], i['fq1'], i['fq2'], f"{wd}/{i['tumor']}", 'run.mrd.sh', argv.nt)
        do_qsub(f"{wd}/{i['tumor']}/run.mrd.sh", argv.nt) if argv.qsub else print(f"write {wd}/{i['tumor']}/run.mrd.sh done.")


