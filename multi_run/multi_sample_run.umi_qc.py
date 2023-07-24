#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
批量刷出多个样本的shell并投递运行
运行umi带qc版本call变异
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


def make_shell(tumor, fq1, fq2, outdir, sh):
    cmd = f"""#!/bin/bash
set -euo pipefail
echo ==== start at `date "+%F  %H:%M:%S"` ====

PL="DNBSEQ"
TUMOR_SM={tumor}
TUMOR_FASTQ_1={fq1}
TUMOR_FASTQ_2={fq2}

# Update with the location of the reference data files
FASTA_DIR="/mnt/share02/zhouyj/database/human_b37_reference"
FASTA="$FASTA_DIR/human_g1k_v37_decoy.chr.fasta"
KNOWN_DBSNP="$FASTA_DIR/dbsnp_138.b37.chr.vcf.gz"
INTERVAL_FILE="/mnt/share03/tumor_pipeline/Somatic/DNA/Pipeline/Somatic_V1/database/bed/depth_bed/NPC69.bed"
TRUTH_VCF="/mnt/share05/clinical_project/projects/blood_tumor/test/zhufu_test/double_UMI/sentieon_qc/truth.normed.vcf.gz"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=/mnt/share02/xuefei/dna/sentieon/sentieon-genomics-202010.04
export SENTIEON_LICENSE=10.90.1.100:8990

#UMI information
READ_STRUCTURE="3M4S+T,3M4S+T"
DUPLEX_UMI="true" #set to "false" if not duplex

# Other settings
NT=$(nproc) #number of threads to use in computation, set to number of cores in the server
START_DIR=$PWD/TNscope_umi

# ******************************************
# Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
cd $WORKDIR

# ******************************************
# 1. Pre-processing of FASTQ containing UMIs
# ******************************************
if [ "$DUPLEX_UMI" = "true" ] ; then
    READ_STRUCTURE="-d $READ_STRUCTURE"
fi
$SENTIEON_INSTALL_DIR/bin/sentieon umi extract $READ_STRUCTURE $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 | \\
  $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C -R "@RG\\tID:$TUMOR_SM\\tSM:$TUMOR_SM\\tPL:$PL" \\
  -t $NT -K 10000000 $FASTA - | tee >($SENTIEON_INSTALL_DIR/bin/sentieon util sort -i - --sam2bam -o $TUMOR_SM.raw.bam) | \\
  $SENTIEON_INSTALL_DIR/bin/sentieon umi consensus -o $TUMOR_SM.umi_consensus.fastq.gz

$SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C -R "@RG\\tID:$TUMOR_SM\\tSM:$TUMOR_SM\\tPL:$PL" \\
  -t $NT -K 10000000 $FASTA $TUMOR_SM.umi_consensus.fastq.gz | \\
  $SENTIEON_INSTALL_DIR/bin/sentieon util sort --umi_post_process --sam2bam -i - -o $TUMOR_SM.umi_consensus.bam

# ******************************************
# 2. Somatic variant calling
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i $TUMOR_SM.umi_consensus.bam \\
    --interval $INTERVAL_FILE \\
    --algo TNscope \\
    --tumor_sample $TUMOR_SM \\
    --dbsnp $KNOWN_DBSNP \\
    --pcr_indel_model NONE \\
    --disable_detector sv \\
    --min_tumor_allele_frac 0.0003 \\
    --filter_t_alt_frac 0.0003 \\
    --max_error_per_read 3 \\
    --min_tumor_lod 3.0 \\
    --min_base_qual 40 \\
    --resample_depth 100000 \\
    --assemble_mode 4 \\
    --given $TRUTH_VCF \\
    $TUMOR_SM.umi_consensus.forcecall.vcf.gz

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i $TUMOR_SM.raw.bam \\
    --interval $INTERVAL_FILE \\
    --algo TNscope \\
    --tumor_sample $TUMOR_SM \\
    --dbsnp $KNOWN_DBSNP \\
    --pcr_indel_model NONE \\
    --disable_detector sv \\
    --min_tumor_allele_frac 0.0003 \\
    --filter_t_alt_frac 0.0003 \\
    --max_error_per_read 3 \\
    --min_tumor_lod 3.0 \\
    --resample_depth 100000 \\
    --assemble_mode 4 \\
    --given $TRUTH_VCF \\
    $TUMOR_SM.raw.forcecall.vcf.gz
    
echo ==== end at `date "+%F  %H:%M:%S"` ====
"""
    with open(f'{outdir}/{sh}', 'w') as w:
        w.write(cmd)


def do_qsub(sh, nt):
    os.chdir(os.path.dirname(sh))
    if subprocess.call(f'qsub -V -cwd -l p={nt} {sh}', shell=True) != 0:  # 运行成功返回0
        print(f'#### {sh} 投递失败！####\\n')


if __name__ == "__main__":
    argv = get_args()
    sample_list = os.path.abspath(argv.sample_list)
    wd = os.path.dirname(sample_list)
    tumors = [t.strip() for t in open(sample_list).readlines()]
    fq1s = [f'{wd}/../naangda/{t}/step0.cutadapt/{t}.clean.R1.fastq.gz' for t in tumors]
    fq2s = [f'{wd}/../naangda/{t}/step0.cutadapt/{t}.clean.R2.fastq.gz' for t in tumors]
    for t, fq1, fq2 in zip(tumors, fq1s, fq2s):
        sh = 'run.umi_qc.sh'
        if not os.path.exists(f'{wd}/{t}'):
            os.mkdir(f'{wd}/{t}')
        make_shell(t, fq1, fq2, f'{wd}/{t}', sh)
        do_qsub(f'{wd}/{t}/{sh}', argv.nt) if argv.qsub else print(f'write {wd}/{t}/{sh} done.')

